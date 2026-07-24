#pragma once

#include "FillBase.hpp"

namespace Slic3r {

class FillArcOverhang final : public Fill
{
public:
    Fill *clone() const override { return new FillArcOverhang(*this); }
    bool  is_self_crossing() override { return false; }
    bool  no_sort() const override { return true; }

protected:
    void _fill_surface_single(const FillParams              &params,
                              unsigned int                   thickness_layers,
                              const std::pair<float, Point> &direction,
                              ExPolygon                      expolygon,
                              Polylines                     &polylines_out) override;

    float _layer_angle(size_t) const override { return 0.f; }
};

} // namespace Slic3r
