
#include <iostream>
#include <string>
#include "operator.h"
#include "contraction.h"
#include "projector.h"
#include "tree.h"
#include "parity.h"
#include "evaluate.h"
#include "metadata.h"
#include "subtraction.h"
#include "amputate.h"
#include "jackknife.h"
#include "msbar.h"
#include "Zq.h"

#include "output.h"

double condition_number(Eigen::MatrixXcd &A) {
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A);
    return svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
}

bool parse_binary_flag(int argc, char **argv, const std::string &flag) {
    for (int i = 1; i < argc; i++) {
        if (argv[i] == flag) {
            return true;
        }
    }
    return false;
}

std::string parse_string_arg(int argc, char **argv, const std::string &flag,
        const std::string &default_value = "") {
    // We use i < argc - 1 so we can access i + 1 safely
    for (int i = 1; i < argc - 1; i++) {
        if (argv[i] == flag) {
            return argv[i + 1];
        }
    }
    return default_value;
}

static void write_Zq_data(const std::string &filename,
        JackknifeDatabase<ComplexD> Zq_ZV,
        JackknifeDatabase<ComplexD> Zq_ZA,
        double ZV, double ZA,
        std::vector<int> trajectories) {
    assert(Zq_ZV.values.size() == Zq_ZA.values.size());
    assert(Zq_ZV.values.size() == trajectories.size());

    std::ofstream file(filename);
    int precision = std::numeric_limits<double>::max_digits10;
    file << std::setprecision(precision);

    file << "ZV " << ZV << std::endl;
    file << "ZA " << ZA << std::endl;

    std::string delim = " ";

    file << "trajectories_excluded" << delim << "Zq_from_ZV"
        << delim << "Zq_from_ZA" << std::endl;

    file << "central_value" << delim << Zq_ZV.central_value()
        << delim << Zq_ZA.central_value() << std::endl;

    for (size_t i = 0; i < trajectories.size(); i++) {
        file << trajectories[i] << delim << Zq_ZV.values[i]
            << delim << Zq_ZA.values[i] << std::endl;
    }
}

void print_usage() {
    std::ostream& os = std::cerr;
    os << "Usage: ./analysis data_directory basis Zq_scheme [options]";
    os << std::endl << std::endl;
    os << "Options:" << std::endl;
    os << "\t--parity <parity>" << std::endl
        << "\t\tChoose parity for operators and projectors (default: even)"
        << " (possible values: even odd)" << std::endl;
    os << "\t--bin-size <num>" << std::endl
        << "\t\tSets the bin size for the jackknife procedure (default: 1)" << std::endl;
    os << "\t--enforce-reality" << std::endl
        << "\t\tEnforce that mixing matrices be real by setting the imaginary part to 0"
        << std::endl;
    os << "\t--reaverage" << std::endl
        << "\t\tReaverage jackknife samples at end of calculation when calculating mean value"
        << std::endl;
    os << "\t--no-zq" << std::endl
        << "\t\tDo not include the factor of Z_q^2 in the Z-factors" << std::endl;
    os << "\t--output <filename>" << std::endl
        << "\t\tFile to output Z-factor jackknife samples" << std::endl;
    os << "\t--zq-output <filename>" << std::endl
        << "\t\tFile to output Zq-factor jackknife samples" << std::endl;

    os << std::endl;
}

void die_usage() {
    print_usage();
    exit(1);
}

int main(int argc, char **argv) {
    if (argc < 4) {
        die_usage();
    }
    std::string dir(argv[1]);
    std::string basis(argv[2]);
    std::string wavefunction(argv[3]);

    std::string parity_str = parse_string_arg(argc, argv, "--parity", "even");
    int bin_size = std::stoi(parse_string_arg(argc, argv, "--bin-size", "1"));

    bool enforce_reality = parse_binary_flag(argc, argv, "--enforce-reality");
    bool reaverage = parse_binary_flag(argc, argv, "--reaverage");
    bool no_zq = parse_binary_flag(argc, argv, "--no-zq");
    std::string output_file = parse_string_arg(argc, argv, "--output");
    std::string Zq_output_file = parse_string_arg(argc, argv, "--zq-output");

    double ZV = std::stod(parse_string_arg(argc, argv, "--known-ZV",  "1.0"));
    double ZA = std::stod(parse_string_arg(argc, argv, "--known-ZA",  "1.0"));

    std::cout << "Input directory: " << dir << std::endl;
    std::cout << "Bin size: " << bin_size << std::endl;

    Parity parity;
    if (parity_str == "even") {
        parity = Parity::EVEN;
    }
    else if (parity_str == "odd") {
        parity = Parity::ODD;
    }
    else {
        std::cerr << "Unkown parity: '" << parity_str << "'" << std::endl;
        std::cerr << "(possibilities are 'even' or 'odd')" << std::endl;
        exit(1);
    }
    std::cout << "Parity: " << parity << std::endl;

    if (enforce_reality) {
        std::cout << "Enforcing reality" << std::endl;
    }
    else {
        std::cout << "Not enforcing reality" << std::endl;
    }

    if (reaverage) {
        std::cout << "Reaveraging for jackknife central value" << std::endl;
    }
    else {
        std::cout << "Not reaveraging for jackknife central value" << std::endl;
    }
    std::cout << "Output file: '" << output_file << "'" << std::endl;
    std::cout << "Zq output file: '" << Zq_output_file << "'" << std::endl;

    OperatorBasis op_basis = OperatorBasis::MASAAKI;
    ProjectorBasis proj_basis = ProjectorBasis::MASAAKI;
    if (basis == "gammamu") {
        std::cout << "Using Greg gamma_mu basis" << std::endl;
        op_basis = OperatorBasis::GREG;
        proj_basis = ProjectorBasis::GREG;
    }
    else if (basis == "qslash") {
        std::cout << "Using Greg qslash basis" << std::endl;
        op_basis = OperatorBasis::GREG;
        proj_basis = ProjectorBasis::GREG_QSLASH;
    }
    
    //Change here
    else if (basis == "Masaaki"){
        std::cout << "Using Masaaki basis" << std::endl;
        op_basis = OperatorBasis::MASAAKI;
        proj_basis = ProjectorBasis::MASAAKI;
    }
    else if (basis == "Masaaki_qslash"){
    	std::cout<<"Using Masaaki qslash"<< std::endl;
    	op_basis = OperatorBasis::MASAAKI_QSLASH;
    	proj_basis = ProjectorBasis::MASAAKI_QSLASH;
    }
    else {
        std::cout << "Unkown operator renormalization scheme '"
            << basis << "'" << std::endl;
        std::cout << "Options are :" << std::endl;
        std::cout<<"'gammamu' for Greg gammamu basis"<<std::endl; 
        std::cout<<"'qslash' for Greg qslash basis"<<std::endl; 
		std::cout<<"'Masaaki' for Masaaki gammamu basis"<<std::endl; 
		std::cout<<"'Masaaki_qslash' for Masaaki qslash basis"<<std::endl; 
        exit(1);
    }
	//Change ends here
	
    std::cout << "Loading data...";
    fflush(stdout);
    RunData run_data = read_data_with_parity(dir, parity);
    std::cout << "Done (Loaded " << run_data.metadata.trajectories.size()
        << " trajectories)." << std::endl;

    std::cout << "Making jackknife database...";
    fflush(stdout);
    JackknifeDatabase<TrajectoryData> jack = make_jackknife(run_data.trajectory_data, bin_size);
    std::cout << "Done." << std::endl;

    std::cout << "Amputating data...";
    fflush(stdout);
    jack.fmap_void([&] (TrajectoryData &data) { amputate_trajectory_data(data); });
    std::cout << "Done." << std::endl;

    std::cout << "Computing Z_q from Z_V" << std::endl;
    std::cout << "Using ZV = " << ZV << std::endl;
    WavefunctionScheme scheme;
    if (wavefunction == "qslash") {
        std::cout << "Using qslash wavefunction renormalization scheme"
            << std::endl;
        scheme = WavefunctionScheme::QSLASH;
    }
    else if (wavefunction == "gammamu") {
        std::cout << "Using gamma_mu wavefunction renormalization scheme"
            << std::endl;
        scheme = WavefunctionScheme::GAMMA_MU;
    }
    else {
        std::cout << "Unkown wavefunction renormalization scheme '"
            << wavefunction << "'" << std::endl;
        std::cout << "Options are 'qslash' and 'gammamu'" << std::endl;
        exit(1);
    }
    JackknifeDatabase<ComplexD> Zq_ZV = jack.fmap([&] (TrajectoryData &data) {
            return Zq_from_ZV(data, run_data.metadata, ZV, scheme);
            });
    std::cout << "Zq = " << Zq_ZV.central_value(reaverage) << " +/- "
        << Zq_ZV.error() << std::endl;

    std::cout << "Computing Z_q from Z_A" << std::endl;
    std::cout << "Using ZA = " << ZA << std::endl;
    JackknifeDatabase<ComplexD> Zq_ZA = jack.fmap([&] (TrajectoryData &data) {
            return Zq_from_ZA(data, run_data.metadata, ZA, scheme);
            });
    std::cout << "Zq = " << Zq_ZA.central_value(reaverage) << " +/- "
        << Zq_ZA.error() << std::endl;

    JackknifeDatabase<ComplexD> Zq_sum = Zq_ZV + Zq_ZA;
    JackknifeDatabase<ComplexD> Zq_avg = Zq_sum.fmap([] (ComplexD z) { return z / 2.0; });
    std::cout << "Average Zq = " << Zq_avg.central_value(reaverage) << " +/- "
        << Zq_avg.error() << std::endl;

    JackknifeDatabase<RealD> Zq = Zq_avg.fmap([] (ComplexD z) {
            return z.real();
            });

    std::cout << "Using Zq = " << Zq.central_value(reaverage) << " +/- "
        << Zq.error() << std::endl;


    JackknifeDatabase<ComplexD> Zq_naive_p2 = jack.fmap([&] (TrajectoryData &data) {
            return Zq_naive_pin(data, run_data.metadata);
            });
    std::cout << "Zq_naive_pin = " << Zq_naive_p2.central_value(reaverage) << " +/- "
        << Zq_naive_p2.error() << std::endl;

    JackknifeDatabase<ComplexD> Zq_naive_p1 = jack.fmap([&] (TrajectoryData &data) {
            return Zq_naive_pout(data, run_data.metadata);
            });
    std::cout << "Zq_naive_pout = " << Zq_naive_p1.central_value(reaverage) << " +/- "
        << Zq_naive_p1.error() << std::endl;

    std::vector<FourQuarkOperator> chiral_basis = get_operators(op_basis);
    std::vector<FourQuarkOperator> external_states = get_external_states(proj_basis);

    std::vector<TwoQuarkOp> subtraction_ops = get_subtraction_operators();
    std::vector<OperatorRepresentation> representations
        = get_representations(proj_basis);

    std::vector<FourQuarkProjector> projectors
        = get_projectors(proj_basis, parity, run_data.metadata);
    std::cout << "Computing tree-level mixings for parity "
        << parity << std::endl;
    Eigen::MatrixXcd tree_mixings = tree_level_mixings(chiral_basis,
            external_states, projectors, representations, parity);
    std::cout << "Done computing tree-level mixings" << std::endl;
    std::cout << "Tree level mixing matrix for parity " << parity << std::endl;
    std::cout << tree_mixings << std::endl;
    std::cout << "(Condition number = " << condition_number(tree_mixings) << ")"
        << std::endl;

    std::vector<TwoQuarkProjector> subtraction_projectors
        = get_subtraction_projectors(parity, run_data.metadata);
    std::vector<OperatorRepresentation> subtraction_reps
        = get_subtraction_reps(proj_basis);
    std::cout << "Computing subtractions...";
    fflush(stdout);
    JackknifeDatabase<Eigen::MatrixXcd>
        subtraction_diff = jack.fmap([&] (TrajectoryData &data) {
                MatrixXcd ret = subtraction_correction(&data,
                        chiral_basis, external_states,
                        projectors, representations,
                        subtraction_ops, subtraction_projectors,
                        subtraction_reps);
                if (enforce_reality)
                    ret = ret.real();
                return ret;
                });
    std::cout << "Done." << std::endl;

    std::cout << "Computing mixings...";
    fflush(stdout);
    JackknifeDatabase<Eigen::MatrixXcd>
        mixings = jack.fmap([&] (TrajectoryData &data) {
                MatrixXcd ret = compute_mixings(data,
                        chiral_basis, external_states, projectors,
                        representations);
                if (enforce_reality)
                    ret = ret.real();
                return ret;
                });
    std::cout << "Done." << std::endl;

    std::cout << "Mixings (real part):" << std::endl;
    std::cout << mixings.central_value(reaverage).real() << std::endl;
    std::cout << "Mixings (imaginary part):" << std::endl;
    std::cout << mixings.central_value(reaverage).imag() << std::endl;

    /*std::cout << "Subtractions (real part):" << std::endl;
    std::cout << subtraction_diff.central_value(reaverage).real() << std::endl;
    std::cout << "Subtractions (imaginary part):" << std::endl;
    std::cout << subtraction_diff.central_value(reaverage).imag() << std::endl;*/

    JackknifeDatabase<Eigen::MatrixXcd> subtracted_mixings = mixings - subtraction_diff;
    std::cout << "Subtracted mixings (real part): " << std::endl;
    std::cout << subtracted_mixings.central_value(reaverage).real() << std::endl;
    std::cout << "Subtracted mixings (imaginary part): " << std::endl;
    std::cout << subtracted_mixings.central_value(reaverage).imag() << std::endl;

    std::cout << "Subtracted mixings Error (real part):" << std::endl;
    std::cout << subtracted_mixings.error().real() << std::endl;
    std::cout << "Subtracted mixings Error (imaginary part):" << std::endl;
    std::cout << subtracted_mixings.error().imag() << std::endl;

    // Lattice -> RI/SMOM conversion factors
    std::cout << "Lat -> RI/SMOM Z factor (No Zq^2):" << std::endl;
    JackknifeDatabase<Eigen::MatrixXcd> zfactor = subtracted_mixings.fmap([&] (Eigen::MatrixXcd &mat) {
            Eigen::MatrixXcd ret = tree_mixings * mat.inverse();
            return ret;
            });
    std::cout << zfactor.central_value(reaverage).real() << std::endl;

    if (!no_zq) {
        std::cout << "Lat -> RI/SMOM Z factor:" << std::endl;
        zfactor *= Zq;
        zfactor *= Zq;
        std::cout << zfactor.central_value(reaverage).real() << std::endl;
    }
    std::cout << "Lat -> RI/SMOM Z factor errors:" << std::endl;
    std::cout << zfactor.error().real() << std::endl;

    // Collect data on trajectory numbers into a vector
    std::vector<int> trajectories;
    trajectories.reserve(run_data.trajectory_data.size());
    for (TrajectoryData &data : run_data.trajectory_data) {
        trajectories.push_back(data.trajectory);
    }


    if (!output_file.empty()) {
        if (op_basis == OperatorBasis::MASAAKI) {
            const int N = 9;
            OutputData<N> out = convert_to_output<N>(zfactor, trajectories);
            out.write_to_file(output_file);
        }
        else {
            const int N = 7;
            OutputData<N> out = convert_to_output<N>(zfactor, trajectories);
            out.write_to_file(output_file);
        }
    }

    if (!Zq_output_file.empty()) {
        write_Zq_data(Zq_output_file, Zq_ZV, Zq_ZA, ZV, ZA, trajectories);
    }

    return 0;
}
