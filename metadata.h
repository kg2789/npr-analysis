#ifndef __METADATA_H__

#define __METADATA_H__

#include <Grid/Grid.h>
#include <Hadrons/Global.hpp>
#include <Hadrons/Application.hpp>

using namespace Grid;
using namespace Hadrons;

class Metadata : Serializable {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string, name,
                                        std::vector<RealD>, p1,
                                        std::vector<RealD>, p2,
                                        std::vector<std::string>, flavours,
                                        std::vector<RealD>, masses,
                                        std::vector<int>, trajectories,
                                        std::vector<int>, lattice_dimensions,
                                        int, num_a2a_vectors);
};

inline void write_metadata(const std::string &filename, const Metadata &meta) {
    XmlWriter writer(filename);
    write(writer, "metadata", meta);
}

inline Metadata read_metadata(const std::string &filename) {
    XmlReader reader(filename);
    Metadata meta;
    read(reader, "metadata", meta);
    return meta;
}

inline double get_a_inv(const Metadata &meta) {
    if (meta.name == "32Ifine") {
        return 3.148;
    }
    else {
        std::cerr << "No value of a_inv found for run '"
            << meta.name << "'" << std::endl;
        exit(1);
    }
    return -1;
}

// Gets the physical scale \mu of the renormalization scale in GeV
inline double get_scale(const Metadata &meta) {
    double a_inv = get_a_inv(meta);
    double scale = 0.0;
    for (int mu = 0; mu < meta.p1.size(); mu++) {
        const int p1mu = meta.p1[mu];
        const int Lmu = meta.lattice_dimensions[mu];
        double p1_physical_mu = 2 * M_PI / Lmu * p1mu;
        scale += p1_physical_mu * p1_physical_mu;
    }
    scale = std::sqrt(scale) * a_inv;
    return scale;
}

#endif
