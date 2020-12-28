//! 3D Checkerboard

#[allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::texture::*;

/// Implements a checkerboard texture via a 3D mapping.
#[derive(Clone)]
pub struct CheckerboardTexture3D<T> {
    /// First texture.
    tex1: ArcTexture<T>,

    /// Second texture.
    tex2: ArcTexture<T>,

    /// 3D mapping.
    mapping: ArcTextureMapping3D,
}

impl<T> CheckerboardTexture3D<T> {
    /// Create a new `CheckerboardTexture3D<T>`.
    ///
    /// * `tex1`      - The first texture.
    /// * `tex2`      - The second texture.
    /// * `mapping`   - The 3D mapping.
    pub fn new(tex1: ArcTexture<T>, tex2: ArcTexture<T>, mapping: ArcTextureMapping3D) -> Self {
        Self {
            tex1: tex1.clone(),
            tex2: tex2.clone(),
            mapping: mapping.clone(),
        }
    }
}

impl<T> Texture<T> for CheckerboardTexture3D<T>
where
    T: Copy,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        // Get the (s, t) mapping for the intersection.
        let TextureMap3DResult { p, .. } = self.mapping.map(si);
        if (p.x.floor() as Int + p.y.floor() as Int + p.z.floor() as Int) % 2 == 0 {
            self.tex1.evaluate(si)
        } else {
            self.tex2.evaluate(si)
        }
    }
}
