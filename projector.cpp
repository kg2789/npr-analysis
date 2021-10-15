#include "projector.h"
#include "metadata.h"

using Grid::Ns;
using Grid::Nc;
using Grid::Nd;

ComplexD FourQuarkProjector::project(SpinColourSpinColourMatrix &op) {
    using Grid::traceIndex;
    SpinColourSpinColourMatrix product = mat * op;
    return TensorRemove(traceIndex<1>(traceIndex<2>(
                traceIndex<3>(traceIndex<4>(product)))));
}

ComplexD TwoQuarkProjector::project(SpinColourMatrix &op) {
    using Grid::traceIndex;
    SpinColourMatrix product = mat * op;
    return TensorRemove(traceIndex<1>(traceIndex<2>(product)));
}

FourQuarkProjector::FourQuarkProjector(const SpinStructure &a,
        const SpinStructure &b,
        ColorStructure structure) {
    mat = vertex_structure(a, b, structure);
}

static std::vector<FourQuarkProjector> get_projectors_greg(Parity parity) {
    std::vector<FourQuarkProjector> projectors;
    projectors.resize(7);

    SpinStructure L = SpinStructure::LEFT;
    SpinStructure R = SpinStructure::RIGHT;
    ColorStructure diag = ColorStructure::COLOR_DIAG;
    ColorStructure mixed = ColorStructure::COLOR_MIXED;

    projectors[0] = FourQuarkProjector(L, L, diag); // P1
    projectors[1] = FourQuarkProjector(L, L, mixed); // P2
    projectors[2] = FourQuarkProjector(L, L, diag); // P3
    projectors[3] = FourQuarkProjector(L, R, diag); // P4
    projectors[4] = FourQuarkProjector(L, R, mixed); // P5
    projectors[5] = FourQuarkProjector(L, R, diag); // P6
    projectors[6] = FourQuarkProjector(L, R, mixed); // P7

    for (FourQuarkProjector &proj: projectors) {
        proj.mat = parity_project(proj.mat, parity);
        // Taken from Greg's code; it appears this minus sign makes the
        // tree-level mixing results the same (otherwise they come out with an
        // overal minus sign)
        if (parity == Parity::ODD) {
            proj.mat = -proj.mat;
        }
    }

    return projectors;
}

// This utility function is taken/translated from Greg's code, NPR_Utils.cpp,
// function BuildQslashProjector
static SpinColourSpinColourMatrix
build_qslash_projector(std::vector<RealD> q, RealD qsq,
        RealD sign, ColorStructure color_structure, Parity parity) {
    assert(sign == +1.0 || sign == -1.0);
    Grid::Gamma g5(Grid::Gamma::Algebra::Gamma5);

    SpinMatrix qslash = vector_slash(q);
    SpinMatrix qslash_g5 = qslash * g5;

    SpinColourSpinColourMatrix even
        = tensor_prod_with_color_structure(qslash, qslash, color_structure)
        + sign * tensor_prod_with_color_structure(qslash_g5, qslash_g5, color_structure);
    SpinColourSpinColourMatrix odd
        = tensor_prod_with_color_structure(qslash_g5, qslash, color_structure)
        + sign * tensor_prod_with_color_structure(qslash, qslash_g5, color_structure);
	
    even = (1.0 / qsq) * even;
    odd = (1.0 / qsq) * odd;

    SpinColourSpinColourMatrix ret;
    switch (parity) {
        case Parity::EVEN:
            ret = even;
            break;
        case Parity::ODD:
            ret = odd;
            break;
        case Parity::BOTH:
            ret = even + odd;
            break;
    }
    return ret;
}

static SpinColourSpinColourMatrix
build_spin_mixed_qslash_projector(std::vector<RealD> q, RealD qsq,
        RealD sign, ColorStructure color_structure, Parity parity) {
    assert(sign == +1.0 || sign == -1.0);
    Grid::Gamma g5(Grid::Gamma::Algebra::Gamma5);

    SpinMatrix qslash = vector_slash(q);
    SpinMatrix qslash_g5 = qslash * g5;

    SpinColourSpinColourMatrix even
        = spin_mixed_tensor_prod_with_color_structure(qslash, qslash, color_structure)
        + sign * spin_mixed_tensor_prod_with_color_structure(qslash_g5, qslash_g5, color_structure);
    SpinColourSpinColourMatrix odd
        = spin_mixed_tensor_prod_with_color_structure(qslash_g5, qslash, color_structure)
        + sign * spin_mixed_tensor_prod_with_color_structure(qslash, qslash_g5, color_structure);
	
    even = (1.0 / qsq) * even;
    odd = (1.0 / qsq) * odd;

    SpinColourSpinColourMatrix ret;
    switch (parity) {
        case Parity::EVEN:
            ret = even;
            break;
        case Parity::ODD:
            ret = odd;
            break;
        case Parity::BOTH:
            ret = even + odd;
            break;
    }
    return ret;
}


// Taken from Greg's code, NPR_Utils.cpp, function
// BuildQslashProjectorSpinColorStructures
static std::vector<FourQuarkProjector>
get_projectors_greg_qslash(Parity parity, Metadata metadata) {
    std::vector<RealD> p1 = metadata.p1;
    std::vector<RealD> p2 = metadata.p2;
    std::vector<RealD> q;
    q.resize(p1.size());
    RealD qsq = 0.0;
    for (int mu = 0; mu < Nd; mu++) {
        int Lmu = metadata.lattice_dimensions[mu];
        p1[mu] *= 2 * M_PI / Lmu;
        p2[mu] *= 2 * M_PI / Lmu;
        q[mu] = p1[mu] - p2[mu];
        qsq += q[mu] * q[mu];
    }

    SpinColourSpinColourMatrix P1plus
        = build_qslash_projector(q, qsq, +1.0,
                ColorStructure::COLOR_DIAG, parity);
    SpinColourSpinColourMatrix P1minus
        = build_qslash_projector(q, qsq, -1.0,
                ColorStructure::COLOR_DIAG, parity);
    SpinColourSpinColourMatrix P2plus
        = build_qslash_projector(q, qsq, +1.0,
                ColorStructure::COLOR_MIXED, parity);
    SpinColourSpinColourMatrix P2minus
        = build_qslash_projector(q, qsq, -1.0,
                ColorStructure::COLOR_MIXED, parity);

    const RealD Nc = 3;

    // All of the citations below are copied directly from Greg's code

    std::vector<SpinColourSpinColourMatrix> ret;
    ret.resize(7);

    // (27,1) projector, Eq. (99) in Lehner & Sturm 1104.4948
    ret[0] = (1.0 / (64 * Nc * (Nc + 1))) * P1plus;

    const double coeff = 1.0 / (32 * Nc * (Nc*Nc - 1));

    // (8, 8) projectors, Eq. 100
    ret[5] = coeff * (Nc * P1minus - P2minus);
    ret[6] = coeff * (-1.0 * P1minus + Nc * P2minus);

    // (8,1) projectors, Eq. 102
    ret[1] = coeff * ((3 * Nc - 2) * P1plus + (2 * Nc - 3) * P2plus);
    ret[2] = coeff * ((2 * Nc - 3) * P1plus + (3 * Nc - 2) * P2plus);
    ret[3] = ret[5];
    ret[4] = ret[6];

    std::vector<FourQuarkProjector> projectors;
    projectors.reserve(ret.size());
    for (SpinColourSpinColourMatrix &mat: ret) {
        projectors.push_back(FourQuarkProjector(mat));
    }

    return projectors;
}

static std::vector<FourQuarkProjector>
get_projectors_daiquian_qslash(Parity parity, Metadata metadata) {
    std::vector<RealD> p1 = metadata.p1;
    std::vector<RealD> p2 = metadata.p2;
    std::vector<RealD> q;
    q.resize(p1.size());
    RealD qsq = 0.0;
    for (int mu = 0; mu < Nd; mu++) {
        int Lmu = metadata.lattice_dimensions[mu];
        p1[mu] *= 2 * M_PI / Lmu;
        p2[mu] *= 2 * M_PI / Lmu;
        q[mu] = p1[mu] - p2[mu];
        qsq += q[mu] * q[mu];
    }

    SpinColourSpinColourMatrix P1plus
        = build_qslash_projector(q, qsq, +1.0,
                ColorStructure::COLOR_DIAG, parity);
    SpinColourSpinColourMatrix P1minus
        = build_qslash_projector(q, qsq, -1.0,
                ColorStructure::COLOR_DIAG, parity);
    SpinColourSpinColourMatrix P2plus
        = build_qslash_projector(q, qsq, +1.0,
                ColorStructure::COLOR_MIXED, parity);
    SpinColourSpinColourMatrix P2minus
        = build_qslash_projector(q, qsq, -1.0,
                ColorStructure::COLOR_MIXED, parity);

    std::vector<SpinColourSpinColourMatrix> ret;
    ret.resize(7);

    // (27, 1)
    ret[0] = P1plus;

    // (8, 1)
    ret[1] = P1plus;
    ret[2] = P2plus;
    ret[3] = P1minus;
    ret[4] = P2minus;

    // (8, 8)
    ret[5] = P1minus;
    ret[6] = P1plus;

    std::vector<FourQuarkProjector> projectors;
    projectors.reserve(ret.size());
    for (SpinColourSpinColourMatrix &mat: ret) {
        projectors.push_back(FourQuarkProjector(mat));
    }

    return projectors;
}

// Extension of Greg's projectors to a 4-quark case
static std::vector<FourQuarkProjector> get_projectors_masaaki(Parity parity) {
    std::vector<FourQuarkProjector> projectors = get_projectors_greg(Parity::BOTH);
    projectors.resize(9);

    SpinStructure L = SpinStructure::LEFT;
    ColorStructure diag = ColorStructure::COLOR_DIAG;
    ColorStructure mixed = ColorStructure::COLOR_MIXED;

    projectors[7] = FourQuarkProjector(L, L, diag);
    projectors[8] = FourQuarkProjector(L, L, mixed);

    for (FourQuarkProjector &proj: projectors) {
        proj.mat = parity_project(proj.mat, parity);
        // Taken from Greg's code; it appears this minus sign makes the
        // tree-level mixing results the same (otherwise they come out with an
        // overal minus sign)
        if (parity == Parity::ODD) {
            proj.mat = -proj.mat;
        }
    }

    return projectors;
}

//Change here
static std::vector<FourQuarkProjector> get_projectors_masaaki_qslash(Parity parity, Metadata metadata) {
	std::vector<RealD> p1 = metadata.p1;
    std::vector<RealD> p2 = metadata.p2;
    std::vector<RealD> q;
    q.resize(p1.size());
    RealD qsq = 0.0;
    for (int mu = 0; mu < Nd; mu++) {
        int Lmu = metadata.lattice_dimensions[mu];
        p1[mu] *= 2 * M_PI / Lmu;
        p2[mu] *= 2 * M_PI / Lmu;
        q[mu] = p1[mu] - p2[mu];
        qsq += q[mu] * q[mu];
    }
    SpinColourSpinColourMatrix P1plus
        = build_qslash_projector(q, qsq, +1.0,
                ColorStructure::COLOR_DIAG, parity);
    SpinColourSpinColourMatrix P1minus
        = build_qslash_projector(q, qsq, -1.0,
                ColorStructure::COLOR_DIAG, parity);
    SpinColourSpinColourMatrix P2plus
        = build_qslash_projector(q, qsq, +1.0,
                ColorStructure::COLOR_MIXED, parity);
    SpinColourSpinColourMatrix P2minus
        = build_qslash_projector(q, qsq, -1.0,
                ColorStructure::COLOR_MIXED, parity);
    
    //I created P3 and P4 for completeness sake
    //I thought I might need them, but never had to
    //However, one can use these as well of course
	SpinColourSpinColourMatrix P3plus
        = build_spin_mixed_qslash_projector(q, qsq, +1.0,
                ColorStructure::COLOR_MIXED, parity);
	SpinColourSpinColourMatrix P3minus
        = build_spin_mixed_qslash_projector(q, qsq, -1.0,
                ColorStructure::COLOR_MIXED, parity);
    SpinColourSpinColourMatrix P4plus
        = build_spin_mixed_qslash_projector(q, qsq, +1.0,
                ColorStructure::COLOR_DIAG, parity);
	SpinColourSpinColourMatrix P4minus
        = build_spin_mixed_qslash_projector(q, qsq, -1.0,
                ColorStructure::COLOR_DIAG, parity);
                                
    //Number of Colors
    const RealD Nc = 3;
    std::vector<SpinColourSpinColourMatrix> ret;
    ret.resize(9);
    
    //See Kshitij notes for all the equations
    
    //There are 2 operators in (84,1) representation, 
    //so should there not be 2 projectors too
    //Calculated in Kshitij notes
    
    //(84,1) projector
    ret[0] = std::sqrt(15)/(80.0 * Nc * (Nc + 1)) * P1plus;
    
    //(84,1) second projector for the two operators (Using same projector)
    ret[1] = ret[0];    
    //ret[1] = P1plus;
	
	//(20,1) projector calculated in Kshitij notes
	//ret[2] = (1.0 / (16 * Nc * (Nc - 1))) * P1plus;
	ret[2] = (1.0/(16.0 * Nc * (Nc - 1))) * P1plus;
	
	ret[3] = (1.0/(16.0 * Nc * (Nc - 1))) * P1plus;
	ret[4] = std::sqrt(3)/(32.0 * (Nc * (Nc + 1))) * P1plus;
	
	ret[5] = (1.0/(16.0 * Nc * Nc)) * P1minus;
	ret[6] = (1.0/(16.0 * Nc * Nc)) * P2minus;
	
	//For Projectors 8 and 9 :
	//Please check Kshitij notes in the individual postings
	//I am missing a factor of 2 in my analytic calculations
	//, but I do not understand where
	//However, this does NOT affect final Z factors
	//so one can ignore this sidenote
	ret[7] = (1.0/(32.0 * Nc * Nc)) * P1minus;
	ret[8] = (1.0/(32.0 * Nc * Nc)) * P2minus;

	
	std::vector<FourQuarkProjector> projectors;
    projectors.reserve(ret.size());
    for (SpinColourSpinColourMatrix &mat: ret) {
        projectors.push_back(FourQuarkProjector(mat));
    }

    return projectors;
	
}
//Change ends here	
	

std::vector<FourQuarkProjector> get_projectors(ProjectorBasis basis, Parity parity,
        Metadata metadata) {
    switch (basis) {
        case ProjectorBasis::GREG: return get_projectors_greg(parity);
        case ProjectorBasis::MASAAKI: return get_projectors_masaaki(parity);
        case ProjectorBasis::GREG_QSLASH: return get_projectors_greg_qslash(parity, metadata);
        //Change here
        case ProjectorBasis::MASAAKI_QSLASH: return get_projectors_masaaki_qslash(parity, metadata);
        //Change ends here
    }
}

static std::vector<TwoQuarkProjector> get_subtraction_projectors(Parity parity,
        std::vector<RealD> &p1, std::vector<RealD> &p2) {
    using Grid::Gamma;
    Gamma g5 = Gamma(Gamma::Algebra::Gamma5);

    std::vector<TwoQuarkProjector> ret;
    ret.resize(3);

    ret[0] = TwoQuarkProjector(spin_colour_identity());
    ret[1] = TwoQuarkProjector(vector_slash_sc(p2));
    ret[2] = TwoQuarkProjector(vector_slash_sc(p1));

    // Project onto parity
    for (int i = 0; i < 3; i++) {
        switch (parity) {
            case Parity::EVEN: break;
            case Parity::ODD:
                ret[i].mat = ret[i].mat * g5;
                break;
            case Parity::BOTH:
                ret[i].mat = ret[i].mat + ret[i].mat * g5;
        }
    }
    return ret;
}

std::vector<TwoQuarkProjector> get_subtraction_projectors(Parity parity,
        Metadata metadata) {
    std::vector<RealD> p1 = metadata.p1;
    std::vector<RealD> p2 = metadata.p2;
    for (int mu = 0; mu < Nd; mu++) {
        int Lmu = metadata.lattice_dimensions[mu];
        p1[mu] *= 2 * M_PI / Lmu;
        p2[mu] *= 2 * M_PI / Lmu;
    }
    return get_subtraction_projectors(parity, p1, p2);
}

// Under some schemes (mainly the qslash scheme) we must enforce a lack of
// mixing between representations as part of the specification of the scheme.
// The reps given here are for both the operators and the projectors, and the
// 'NONE' representation is a stand-in for when we don't care about
// representations. As an example, the scheme given for GREG_QSLASH will force
// the mixing matrix to look like this:
//
// x 0 0 0 0 0 0
// 0 x x x x 0 0
// 0 x x x x 0 0
// 0 x x x x 0 0
// 0 x x x x 0 0
// 0 0 0 0 0 x x
// 0 0 0 0 0 x x
//
// Where 'x' is some value. Meanwhile for the basis GREG every entry of the
// matrix could be nonzero.
//
// In the future we might need to specify representations separately for
// operators and projectors; however, every scheme we have uses the same set of
// representations for both operators and projectors, so there's no reason to
// make this distinction.
std::vector<OperatorRepresentation> get_representations(ProjectorBasis basis) {
    switch (basis) {
        case ProjectorBasis::GREG:
            return std::vector<OperatorRepresentation>(7,
                    OperatorRepresentation::NONE);
        case ProjectorBasis::MASAAKI:
            return std::vector<OperatorRepresentation>(9,
                    OperatorRepresentation::NONE);
        case ProjectorBasis::GREG_QSLASH:
            return {
                OperatorRepresentation::REP27_1,
                OperatorRepresentation::REP8_1,
                OperatorRepresentation::REP8_1,
                OperatorRepresentation::REP8_1,
                OperatorRepresentation::REP8_1,
                OperatorRepresentation::REP8_8,
                OperatorRepresentation::REP8_8
            };
        //Change here
        case ProjectorBasis::MASAAKI_QSLASH:
        	return{
        		OperatorRepresentation::REP84_1,
        		OperatorRepresentation::REP84_1,
        		OperatorRepresentation::REP20_1,
        		OperatorRepresentation::REP15_1,
        		OperatorRepresentation::REP15_1,
        		OperatorRepresentation::REP15_1,
				OperatorRepresentation::REP15_1,
        		OperatorRepresentation::REP15_15,
        		OperatorRepresentation::REP15_15
        	};
		//Change ends here
        default:
            std::cerr << "Error: unknown projector basis" << std::endl;
            exit(1);
    }
}
