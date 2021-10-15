#include "io.h"

void write_trajectory(const std::string& dir, TrajectoryData* data,
        Metadata &metadata, int trajectory, Format format) {
    switch (format) {
        case Format::RABBOTT:
            return write_trajectory_rabbott(dir, data, trajectory);
        case Format::GREG:
            return write_trajectory_greg(dir, data, metadata, trajectory);
        default:
            std::cerr << "Error: format " << format << " not supported for writing";
            exit(1);
    }
}

void read_trajectory(const std::string& dir, TrajectoryData *data,
        Metadata &metadata, int trajectory, Format format) {
    data->trajectory = trajectory;
    switch (format) {
        case Format::RAW:
            return read_trajectory_raw(dir, data, trajectory);
        case Format::RAW_V2:
            return read_trajectory_raw_v2(dir, data, metadata, trajectory);
        case Format::RAW_WITH_MIXED:
            return read_trajectory_raw_with_mixed(dir, data, trajectory);
        case Format::RABBOTT:
            return read_trajectory_rabbott(dir, data, trajectory);
        case Format::GREG:
            return read_trajectory_greg(dir, data, metadata, trajectory);
        case Format::ZIYUAN:
            return read_trajectory_ziyuan(dir, data, metadata, trajectory);
        default:
            std::cerr << "Error: format " << format
                << " not supported for reading";
            exit(1);
    }
}

