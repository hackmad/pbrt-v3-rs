//! GUI for displaying the rendered result

use std::sync::{OnceLock, RwLock};

use pixels::{Pixels, SurfaceTexture};
use tao::{
    dpi::LogicalSize,
    event::{ElementState, Event, KeyEvent, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    keyboard::Key,
    window::{Window, WindowBuilder},
};

/// The GUI.
pub struct Gui {
    /// Event loop.
    event_loop: EventLoop<()>,
}

/// The window inner size in pixels.
pub const WINDOW_WIDTH: u32 = 400;
pub const WINDOW_HEIGHT: u32 = 400;

/// The global window.
pub static WINDOW: OnceLock<Window> = OnceLock::new();

/// The global window pixel frame buffer.
pub static WINDOW_PIXELS: OnceLock<RwLock<Pixels>> = OnceLock::new();

impl Gui {
    /// Build a new `Gui`.
    pub fn build() -> Result<Self, String> {
        // Create a new event loop for the application.
        let event_loop = EventLoop::new();

        // Create a new window.
        let win = WINDOW.get_or_init(|| {
            WindowBuilder::new()
                .with_title("PBRT v3 (Rust)")
                .with_inner_size(LogicalSize::new(WINDOW_WIDTH, WINDOW_HEIGHT))
                .with_resizable(true)
                .build(&event_loop)
                .expect("Unable to create window")
        });

        let inner_size = win.inner_size();

        // Create a surface texture that uses the logical inner size to render to the entire window's inner dimensions.
        let surface_texture = SurfaceTexture::new(inner_size.width, inner_size.height, win);

        // Create pixel frame buffer that matches rendered image dimensions that will be used to display it in the
        // window.
        WINDOW_PIXELS.get_or_init(|| {
            RwLock::new(
                Pixels::new(WINDOW_WIDTH, WINDOW_HEIGHT, surface_texture)
                    .expect("Unable to create pixel frame buffer for window"),
            )
        });

        Ok(Self { event_loop })
    }

    /// Run the event loop displaying the GUI window until it is closed or some error occurs.
    ///
    /// NOTE: This consumes `self` to avoid life time issues with references and also it runs indefinitely so it can
    /// only ever be called once.
    pub fn run(self) -> ! {
        let Self { event_loop } = self;

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
                    _ => (),
                },
                Event::RedrawRequested(_) => {
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
                _ => (),
            }
        })
    }
}
