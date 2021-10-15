#ifndef __IO_H__

#define __IO_H__

#include <sys/stat.h>

#include "format.h"
#include "../data.h"

/* All read/write functions return success of the operation */
void write_trajectory(const std::string& dir, TrajectoryData* data,
        Metadata &metadata, int trajectory, Format format);
void read_trajectory(const std::string& dir, TrajectoryData* data,
      Metadata &metadata, int trajectory, Format format);

void write_trajectory_rabbott(const std::string& dir, TrajectoryData* data, int trajectory);
void read_trajectory_rabbott(const std::string& dir, TrajectoryData* data, int trajectory);

void read_trajectory_raw(const std::string& dir, TrajectoryData* data, int trajectory);
void read_trajectory_raw_v2(const std::string& dir, TrajectoryData* data, Metadata &metadata, int trajectory);
void read_trajectory_raw_with_mixed(const std::string& dir, TrajectoryData* data, int trajectory);

void write_trajectory_greg(const std::string& dir, TrajectoryData* data,
        Metadata &metadata, int trajectory);
void read_trajectory_greg(const std::string& dir, TrajectoryData* data,
        Metadata &metadata, int trajectory);

void read_trajectory_ziyuan(const std::string& dir, TrajectoryData* data,
        Metadata &metadata, int trajectory);

static inline void make_dir(const std::string& dirname) {
    int res = mkdir(dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (res != 0 && errno != EEXIST) {
        std::string s = "Error making directory '" + dirname + "'";
        perror(s.c_str());
        exit(1);
    }
}
#endif
