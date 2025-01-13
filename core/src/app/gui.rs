//! Window for displaying the rendered result

use std::{sync::{Arc, OnceLock}, thread, time::Duration};

use pixels::{Pixels, SurfaceTexture};
use winit::{
    application::ApplicationHandler,
    dpi::{LogicalSize, PhysicalSize},
    error::EventLoopError,
    event::{ElementState, KeyEvent, WindowEvent},
    event_loop::{ActiveEventLoop, EventLoop, EventLoopProxy},
    keyboard::{Key, NamedKey},
    window::Window,
};

use crate::geometry::Bounds2i;

use super::OPTIONS;

/// The starting window inner size in pixels.
const DEFAULT_WINDOW_WIDTH: u32 = 400;
const DEFAULT_WINDOW_HEIGHT: u32 = 400;

/// User events for the render loop.
#[derive(Debug, Clone, PartialEq)]
enum UserEvent{
    // Render the preview.
    RenderPreview {
        /// The tile pixel data.
        pixels: Vec<u8>,

        /// The starting offset of the merged pixel in the film.
        starting_merge_pixel: usize,

        /// Cropped width of the subset of the image to render.
        cropped_pixel_width: usize,

        /// The tile pixel bounds.
        tile_pixel_bounds: Bounds2i,

        /// The tile pixel width.
        tile_width: usize,
    },
    // Clear the preview window.
    ClearPreview,

    /// Resize the preview window.
    ResizePreview {
        /// The new image size for the preview.
        new_pixel_size: LogicalSize<u32>, 
    },
}

/// This proxy will be used to trigger custom events from the render loop to the winit application window.
static EVENT_LOOP_PROXY: OnceLock<EventLoopProxy<UserEvent>> = OnceLock::new();

/// The winit application.
struct App {
    /// The preview window.
    window: Option<Arc<Window>>,

    /// The preview image pixels.
    pixels: Option<Pixels<'static>>,

    /// The preview image pixel dimensions.
    pixel_size: LogicalSize<u32>,

    /// The inner dimensions of the preview window.
    window_inner_size: PhysicalSize<u32>,
}

impl App {
    /// Render the preview image to the window.
    fn render(&self) -> Result<(), String> {
        self.pixels.as_ref().map_or(Ok(()), |pixels| pixels.render())
            .map_err(|err| format!("{}", err))
    }

    /// Resize the preview image.
    ///
    /// * `pixel_size`        - The dimensions of the preview image.
    /// * `window_inner_size` - The inner dimensions of the preview window.
    fn resize_pixels(
        &mut self,
        pixel_size: LogicalSize<u32>,
        window_inner_size: PhysicalSize<u32>,
    ) -> Result<(), String> {
        // Render only if the application has initialized and we have pixels and window.
        self.pixels.as_mut().map_or(Ok(()), |pixels| {
            // Resize the pixel surface texture to fit the windows inner dimensions.
            match pixels.resize_surface(window_inner_size.width, window_inner_size.height) {
                Ok(()) => {
                    // Resize the pixel image buffer.
                    match pixels.resize_buffer(pixel_size.width, pixel_size.height) {
                        Ok(()) => {
                            // Store the new sizes.
                            self.pixel_size = pixel_size;
                            self.window_inner_size = window_inner_size;

                            // Request a redraw.
                            self.window.as_ref().map(|window| window.request_redraw());
                            Ok(())
                        }
                        Err(err) => Err(format!("pixels.resize_buffer() failed.\n{}", err)),
                    }
                }
                Err(err) => Err(format!("pixels.resize_surface() failed to resize frame buffer surface.\n{}", err)),
            }
        })
    }
}

impl Default for App {
    /// Returns the "default value" for `App` initialized to the default dimensions.
    fn default() -> Self {
        Self {
            window: None,
            pixels: None,
            pixel_size: LogicalSize::new(DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT),
            window_inner_size: PhysicalSize::new(DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT),
        }
    }
}

impl ApplicationHandler<UserEvent> for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        // Create a new window.
        let window_attributes = Window::default_attributes()
            .with_title("PBRT v3 (Rust)")
            .with_inner_size(self.window_inner_size)
            .with_resizable(true);

        let window = Arc::new(event_loop.create_window(window_attributes).expect("Unable to create window"));

        // Save the inner dimensions of the preview window.
        let window_inner_size = window.inner_size();

        // Create a surface texture that uses the logical inner size to render to the entire window's inner
        // dimensions.
        let surface_texture = SurfaceTexture::new(
            window_inner_size.width,
            window_inner_size.height,
            Arc::clone(&window),
        );

        // Create pixel frame buffer that matches rendered image dimensions that will be used to display it
        // in the window.
        let pixels = Pixels::new(self.pixel_size.width, self.pixel_size.height, surface_texture)
            .expect("Unable to create pixel frame buffer for window");

        self.window = Some(Arc::clone(&window));
        self.pixels = Some(pixels);
        self.window_inner_size = window_inner_size;
    }

    fn window_event(
        &mut self,
        event_loop: &ActiveEventLoop,
        _window_id: winit::window::WindowId,
        event: winit::event::WindowEvent,
    ) {
        match event {
            WindowEvent::CloseRequested => {
                println!("The close button was pressed; stopping");
                event_loop.exit();
            }

            WindowEvent::RedrawRequested => {
                match self.render() {
                    Ok(()) => (),
                    Err(err) => {
                        eprintln!("Error redrawing pixels {}", err);
                        event_loop.exit();
                    }
                }
                self.window.as_ref().map(|window| window.request_redraw());
            }
            
            WindowEvent::Resized(new_window_inner_size) => {
                match self.resize_pixels(self.pixel_size, new_window_inner_size) {
                    Ok(()) => (),
                    Err(err) => {
                        eprintln!("Error resizing window {}", err);
                        event_loop.exit();
                    }
                }
            }

            WindowEvent::KeyboardInput {
                event:
                    KeyEvent {
                        logical_key: key,
                        state: ElementState::Pressed,
                        ..
                    },
                ..
            } => match key {
                Key::Named(NamedKey::Escape) => {
                    println!("Escape key was pressed; stopping");
                    event_loop.exit();
                }
                _ => (),
            },

            _ => (),
        }
    }

    fn user_event(&mut self, event_loop: &ActiveEventLoop, event: UserEvent) { 
        match event {
            UserEvent::RenderPreview {
                pixels,
                starting_merge_pixel,
                cropped_pixel_width,
                tile_pixel_bounds,
                tile_width,
            } => {
                self.pixels.as_mut().map(|window_pixels| {
                    let frame = window_pixels.frame_mut();

                    // Some of this logic is duplicated from Film::merge_film_tile().
                    let mut x = tile_pixel_bounds.p_min.x;
                    let mut merge_pixel = starting_merge_pixel;

                    for i in (0..pixels.len()).step_by(4) { // RGBA
                        // Account for window dimensions not matching rendered image dimensions by scaling the
                        // tile appropriately.
                        let window_pixel_x = merge_pixel % cropped_pixel_width as usize;
                        let window_pixel_y = merge_pixel / cropped_pixel_width as usize;
                        let window_pixel_offset = (window_pixel_y * cropped_pixel_width as usize + window_pixel_x) * 4;

                        frame[window_pixel_offset + 0] = pixels[i + 0];
                        frame[window_pixel_offset + 1] = pixels[i + 1];
                        frame[window_pixel_offset + 2] = pixels[i + 2];
                        frame[window_pixel_offset + 3] = pixels[i + 3];

                        x += 1;
                        merge_pixel += 1;

                        // We need to jump to the next line if we've exceeded the horizontal tile bounds.
                        if x == tile_pixel_bounds.p_max.x {
                            x = tile_pixel_bounds.p_min.x;
                            merge_pixel += cropped_pixel_width - tile_width;
                        }
                    }

                    self.window.as_ref().map(|window| window.request_redraw());
                });
            }

            UserEvent::ClearPreview => {
                if let Some(window_pixels) = self.pixels.as_mut() {
                    // Clear the preview image to black.
                    for pixel in window_pixels.frame_mut().chunks_exact_mut(4) {
                        pixel[0] = 0x00; // R
                        pixel[1] = 0x00; // G
                        pixel[2] = 0x00; // B
                        pixel[3] = 0xff; // A
                    }

                    self.window.as_ref().map(|window| window.request_redraw());
                }
            }

            UserEvent::ResizePreview { new_pixel_size } => {
                let resize_status = self.window.as_ref()
                    .map(|window| window.request_inner_size(new_pixel_size).is_some());

                match resize_status {
                    Some(true) => {
                        // If we have a window retrieve its inner size and resize the pixels.
                        if let Some(window_inner_size) = self.window.as_ref().map(|window| window.inner_size()) {
                            if let Err(err) = self.resize_pixels(new_pixel_size, window_inner_size) {
                                eprintln!("{}", err);
                                event_loop.exit()
                            }
                        } 
                    }
                    Some(false) => {
                        // Wait for next WindowEvent::Resized. Store new pixel_size for now.
                        self.pixel_size = new_pixel_size;
                    }
                    None => (), // No window yet
                }
            }
        }
    }
}

/// Run the event loop displaying a window until it is closed or some error occurs.
pub fn run_event_loop() -> Result<(), EventLoopError> {
    eprintln!("Creating event loop");
    let event_loop = EventLoop::<UserEvent>::with_user_event().build().expect("Unable to create event loop");

    eprintln!("Creating up event loop proxy");
    EVENT_LOOP_PROXY.get_or_init(|| event_loop.create_proxy());

    eprintln!("Running winit app");
    let mut app = App::default();
    event_loop.run_app(&mut app)
}

/// Send a user event to the event loop.
///
/// * `event` - The user event to send.
fn send_user_event(event: UserEvent) {
    if OPTIONS.show_gui {
        // The rendering is done a different thread. We could end up here before the event loop is created. So just 
        // check and wait until event loop is ready. This loop will execute only once when the first scene starts 
        // processing.
        while EVENT_LOOP_PROXY.get().is_none() {
            thread::sleep(Duration::from_millis(100));
        }
        EVENT_LOOP_PROXY.get().map(|proxy| proxy.send_event(event));
    }
}

/// Set the window size. This function will wait for window and its frame buffer to be initialized.
///
/// * `bounds` - Crop window of the subset of the image to render.
pub fn set_preview_window_size(bounds: Bounds2i) {
    let size = bounds.diagonal();
    send_user_event(UserEvent::ResizePreview {
        new_pixel_size: LogicalSize::new(size.x as u32, size.y as u32),
    });
}

/// Clear the preview window to black.
pub fn clear_preview_window() {
    send_user_event(UserEvent::ClearPreview);
}

/// Render the preview window's pixel frame buffer.
///
/// * `pixels`               - The tile to render.
/// * `starting_merge_pixel` - The starting offset of the merged pixel in the film.
/// * `cropped_pixel_width`  - Cropped width of the subset of the image to render.
/// * `tile_pixel_bounds`    - The tile pixel bounds.
pub fn render_preview(
    pixels: &[u8],
    starting_merge_pixel: usize,
    cropped_pixel_width: usize,
    tile_pixel_bounds: Bounds2i,
    tile_width: usize,
) {
    send_user_event(UserEvent::RenderPreview {
        pixels: Vec::from(pixels), 
        starting_merge_pixel,
        cropped_pixel_width,
        tile_pixel_bounds,
        tile_width,
    });
}
