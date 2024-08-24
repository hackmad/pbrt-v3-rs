//! Window for displaying the rendered result

use std::{
    sync::{OnceLock, RwLock},
    thread,
    time::Duration,
};

use pixels::{Pixels, SurfaceTexture};
use tao::{
    dpi::{LogicalSize, PhysicalSize},
    event::{ElementState, Event, KeyEvent, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    keyboard::Key,
    window::{Window, WindowBuilder},
};

use crate::{app::OPTIONS, geometry::Bounds2i};

/// The starting window inner size in pixels.
const WINDOW_WIDTH: u32 = 400;
const WINDOW_HEIGHT: u32 = 400;

/// The preview window.
static WINDOW: OnceLock<Window> = OnceLock::new();

/// The pixel frame buffer of the preview window.
static WINDOW_PIXELS: OnceLock<RwLock<Pixels>> = OnceLock::new();

/// Set the window size. This function will wait for window and its frame buffer to be initialized.
///
/// * `bounds` - Crop window of the subset of the image to render.
pub fn set_preview_window_size(bounds: Bounds2i) {
    if OPTIONS.show_gui {
        let size = bounds.diagonal();
        let width = size.x as u32;
        let height = size.y as u32;

        eprintln!("Setting up preview window {}x{}", width, height);

        loop {
            match (
                WINDOW.get(),
                WINDOW_PIXELS.get().map(|pixels| pixels.write().ok()).flatten(),
            ) {
                (Some(window), Some(mut pixels)) => {
                    let inner_size = LogicalSize::new(width, height);
                    window.set_inner_size(inner_size);

                    pixels.resize_surface(inner_size.width, inner_size.height).unwrap();
                    pixels.resize_buffer(width, height).unwrap();

                    break;
                }
                _ => {
                    eprintln!("\rWaiting for window creation");
                    thread::sleep(Duration::from_secs(1));
                }
            }
        }
    }
}

/// Clear the preview window to black.
pub fn clear_preview_window() {
    if OPTIONS.show_gui {
        WINDOW_PIXELS
            .get()
            .map(|pixels| pixels.write().ok())
            .flatten()
            .map(|mut pixels| {
                for pixel in pixels.frame_mut().chunks_exact_mut(4) {
                    pixel[0] = 0x00; // R
                    pixel[1] = 0x00; // G
                    pixel[2] = 0x00; // B
                    pixel[3] = 0xff; // A
                }
            });
    }
}

/// Render the preview window's pixel frame buffer.
///
/// * `f` - Callback that will receive the RGBA byte values of the frame buffer.
pub fn render_preview<F>(f: F)
where
    F: Fn(&mut [u8]),
{
    if OPTIONS.show_gui {
        WINDOW_PIXELS
            .get()
            .map(|wp| wp.write().ok())
            .flatten()
            .map(|mut window_pixels| f(window_pixels.frame_mut()));
    }
}

/// Send request to preview window's event loop to redraw it.
pub fn request_preview_redraw() {
    WINDOW.get().map(|w| w.request_redraw());
}

/// Run the event loop displaying a window until it is closed or some error occurs.
///
/// TODO - When the event loop terminates find a way to signal threads to terminate and stop the rendering process.
pub fn run_event_loop() -> ! {
    // Create a new event loop for the application.
    let event_loop = EventLoop::new();

    // Create a new window with a starting size.
    let window = WINDOW.get_or_init(|| {
        WindowBuilder::new()
            .with_title("PBRT v3 (Rust)")
            .with_inner_size(LogicalSize::new(WINDOW_WIDTH, WINDOW_HEIGHT))
            .with_resizable(true)
            .build(&event_loop)
            .expect("Unable to create window")
    });

    let inner_size = window.inner_size();

    // Create a surface texture that uses the logical inner size to render to the entire window's inner dimensions.
    let surface_texture = SurfaceTexture::new(inner_size.width, inner_size.height, window);

    // Create pixel frame buffer that matches rendered image dimensions that will be used to display it in the
    // window.
    WINDOW_PIXELS.get_or_init(|| {
        RwLock::new(
            Pixels::new(WINDOW_WIDTH, WINDOW_HEIGHT, surface_texture)
                .expect("Unable to create pixel frame buffer for window"),
        )
    });

    event_loop.run(move |event, _, control_flow| {
        //println!("{:?}", event);
        *control_flow = ControlFlow::Wait;

        match event {
            Event::WindowEvent { event, .. } => match event {
                // When window is closed or destroyed or Escape key is pressed, stop rendering.
                WindowEvent::CloseRequested
                | WindowEvent::Destroyed
                | WindowEvent::KeyboardInput {
                    event:
                        KeyEvent {
                            logical_key: Key::Escape,
                            state: ElementState::Released,
                            ..
                        },
                    ..
                } => {
                    eprintln!("Exiting application.");
                    *control_flow = ControlFlow::Exit;
                }
                WindowEvent::Resized(size) => resize(size, control_flow),
                _ => (),
            },
            Event::RedrawRequested(_) => redraw(control_flow),
            _ => (),
        }
    })
}

/// Called when the window is resized. It will resize the pixel surface to ensure it can scale to fit. Note that pixels
/// library does its best to resize so its not perfect.
///
/// * `size`         - The physical size of the window.
/// * `control_flow` - Used to terminate the event loop in case of errors.
fn resize(size: PhysicalSize<u32>, control_flow: &mut ControlFlow) {
    match WINDOW_PIXELS.get().map(|pixels| pixels.write()) {
        Some(Ok(mut pixels)) => match pixels.resize_surface(size.width, size.height) {
            Ok(()) => {
                // Use redraw event handler to update window. Don't call `redraw()` because it will lock up the app.
                WINDOW.get().map(|window| window.request_redraw());
            }
            Err(err) => {
                println!("pixels.render() failed to resize pixel frame buffer surface.\n{}", err);
                *control_flow = ControlFlow::Exit;
            }
        },
        Some(Err(err)) => {
            println!("pixels.render() failed with error.\n{}", err);
        }
        None => {
            println!("pixels.render() failed. Unable to get pixel frame buffer");
        }
    }
}

/// Called when the window needs to be redrawn.
///
/// * `control_flow` - Used to terminate the event loop in case of errors.
fn redraw(control_flow: &mut ControlFlow) {
    // Draw the pixel frame buffer to the window. If there are errors show the error and stop rendering.
    match WINDOW_PIXELS.get().map(|pixels| pixels.read()) {
        Some(Ok(pixels)) => match pixels.render() {
            Ok(()) => {}
            Err(err) => {
                println!("pixels.render() failed with error.\n{}", err);
                *control_flow = ControlFlow::Exit;
            }
        },
        Some(Err(err)) => {
            println!("pixels.render() failed with error.\n{}", err);
            *control_flow = ControlFlow::Exit;
        }
        None => {
            println!("pixels.render() failed. Unable to get pixel frame buffer");
            *control_flow = ControlFlow::Exit;
        }
    }
}
