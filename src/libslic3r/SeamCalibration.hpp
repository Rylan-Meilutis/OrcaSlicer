#pragma once

#include <vector>

namespace Slic3r {
class Model;
class DynamicPrintConfig;

enum class SeamCalibrationKind { Length, StartHeight, Flow, Speed, InnerEntry };
struct SeamCalibrationParams {
    SeamCalibrationKind kind = SeamCalibrationKind::Length;
    double start = 5., step = 5.;
    int count = 4;
    double nozzle = 0.4;
};

// Eight-mm bands: a conventional seam control, then one swept setting per band.
// InnerEntry instead creates matched conventional / inner-entry coupons.
void make_seam_calibration(Model &model, const DynamicPrintConfig &process, const SeamCalibrationParams &params);
std::vector<double> seam_calibration_values(const SeamCalibrationParams &params);
}
