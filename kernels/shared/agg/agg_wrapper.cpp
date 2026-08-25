#include "agg_path_storage.h"
#include "agg_pixfmt_gray.h"
#include "agg_rasterizer_scanline_aa.h"
#include "agg_renderer_base.h"
#include "agg_rendering_buffer.h"
#include "agg_scanline_u.h"
#include "agg_renderer_scanline.h"


extern "C" {


// Draw anti-aliased polygons into the provided buffer.
// The input buffer is assumed to be row-major.
void agg_aa(
    unsigned char* buf,
    int nrows,
    int ncols,
    const double* vx,
    const double* vy,
    int n
) {
    if (n < 3) return;

    agg::path_storage path;
    path.move_to(vx[0], vy[0]);
    for (int i = 1; i < n; ++i) {
        path.line_to(vx[i], vy[i]);
    }
    path.close_polygon();

    agg::rendering_buffer rbuf(
        buf,
        static_cast<unsigned>(ncols),
        static_cast<unsigned>(nrows),
        static_cast<int>(ncols)
    );

    agg::pixfmt_gray8 pixfmt(rbuf);
    agg::renderer_base<decltype(pixfmt)> ren_base(pixfmt);
    ren_base.clear(agg::gray8(0, 255));

    agg::rasterizer_scanline_aa<> ras;
    agg::scanline_u8 sl;
    ras.add_path(path);
    agg::render_scanlines_aa_solid(ras, sl, ren_base, agg::gray8(255, 255));
}


} // extern "C"
