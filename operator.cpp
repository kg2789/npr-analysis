#include "operator.h"
#include "evaluate.h"

FourQuarkOperator tensor_prod(const BilinearOperator &b1, const BilinearOperator &b2, const ColorStructure structure)
{
    FourQuarkOperator ret;
    for (auto &i:  b1.term_coefficients) {
        const BilinearTerm &t1 = i.first;
        const RealD &coeff1 = i.second;
        for (auto &j: b2.term_coefficients) {
            const BilinearTerm &t2 = j.first;
            const RealD &coeff2 = j.second;

            FourQuarkTerm new_term = tensor_prod(t1, t2, structure);
            ret += coeff1 * coeff2 * FourQuarkOperator(new_term);
        }
    }
    return ret;
}

// Only makes 3-quark operators
static std::vector<FourQuarkOperator> get_historical_operators() {
    std::vector<FourQuarkOperator> Q;
    Q.resize(10);
    
    BilinearTerm sd_left((Quark::SBAR), Quark::D, SpinStructure::LEFT);
    BilinearTerm uu_left((Quark::UBAR), Quark::U, SpinStructure::LEFT);

    Q[0] = tensor_prod(sd_left, uu_left, ColorStructure::COLOR_DIAG);   // Q1
    Q[1] = tensor_prod(sd_left, uu_left, ColorStructure::COLOR_MIXED);  // Q2

    BilinearOperator qq_left;
    BilinearOperator qq_right;
    std::vector<Quark> quarks = {{ Quark::U, Quark::D, Quark::S }};
    for(Quark &q: quarks) {
        qq_left += BilinearTerm(antiparticle(q), q, SpinStructure::LEFT);
        qq_right += BilinearTerm(antiparticle(q), q, SpinStructure::RIGHT);
    }

    Q[2] = tensor_prod(sd_left, qq_left, ColorStructure::COLOR_DIAG);  // Q3
    Q[3] = tensor_prod(sd_left, qq_left, ColorStructure::COLOR_MIXED); // Q4
    Q[4] = tensor_prod(sd_left, qq_right, ColorStructure::COLOR_DIAG);  // Q5
    Q[5] = tensor_prod(sd_left, qq_right, ColorStructure::COLOR_MIXED); // Q6

    std::vector<RealD> charges = {{ 2.0 / 3, -1.0 / 3, -1.0 / 3 }};
    BilinearOperator qq_charged_left;
    BilinearOperator qq_charged_right;
    for (int i = 0; i < quarks.size(); i++) {
        Quark &q = quarks[i];
        RealD eq = charges[i];

        qq_charged_left += eq * BilinearOperator(BilinearTerm(antiparticle(q), q, SpinStructure::LEFT));
        qq_charged_right += eq * BilinearOperator(BilinearTerm(antiparticle(q), q, SpinStructure::RIGHT));
    }

    Q[6] = (3.0 / 2) * tensor_prod(sd_left, qq_charged_right, ColorStructure::COLOR_DIAG); // Q7
    Q[7] = (3.0 / 2) * tensor_prod(sd_left, qq_charged_right, ColorStructure::COLOR_MIXED); // Q8
    Q[8] = (3.0 / 2) * tensor_prod(sd_left, qq_charged_left, ColorStructure::COLOR_DIAG); // Q9
    Q[9] = (3.0 / 2) * tensor_prod(sd_left, qq_charged_left, ColorStructure::COLOR_MIXED); // Q10

    return Q;
}

static std::vector<FourQuarkOperator> get_operators_greg() {
    std::vector<FourQuarkOperator> Q = get_historical_operators();
    std::vector<FourQuarkOperator> Qpp;
    Qpp.resize(7);
    FourQuarkOperator Q1 = Q[0];
    FourQuarkOperator Q2 = Q[1];
    FourQuarkOperator Q3 = Q[2];
    FourQuarkOperator Q4 = Q[3];
    FourQuarkOperator Q5 = Q[4];
    FourQuarkOperator Q6 = Q[5];
    FourQuarkOperator Q7 = Q[6];
    FourQuarkOperator Q8 = Q[7];

    Qpp[0] = 3 * Q1 + 2 * Q2 - Q3;
    Qpp[1] = (1.0 / 5) * (2 * Q1 - 2 * Q2 + Q3);
    Qpp[2] = (1.0 / 5) * (-3 * Q1 + 3 * Q2 + Q3);
    Qpp[3] = Q5;
    Qpp[4] = Q6;
    Qpp[5] = Q7;
    Qpp[6] = Q8;

    return Qpp;
}

//copied from Ryans get_operators_masaaki, where I only return first 3 operators since I am testing those
//All comments are same from Ryan's corresponding function
static std::vector<FourQuarkOperator> get_operators_masaaki_qslash(){
	std::vector<FourQuarkOperator> Ppp;
	Ppp.resize(9);
	BilinearTerm sd_left((Quark::SBAR), Quark::D, SpinStructure::LEFT);
    BilinearTerm uu_left((Quark::UBAR), Quark::U, SpinStructure::LEFT);
    BilinearTerm dd_left((Quark::DBAR), Quark::D, SpinStructure::LEFT);
    BilinearTerm ss_left((Quark::SBAR), Quark::S, SpinStructure::LEFT);
    BilinearTerm cc_left((Quark::CBAR), Quark::C, SpinStructure::LEFT);
    BilinearTerm uu_right((Quark::UBAR), Quark::U, SpinStructure::RIGHT);
    BilinearTerm dd_right((Quark::DBAR), Quark::D, SpinStructure::RIGHT);
    BilinearTerm ss_right((Quark::SBAR), Quark::S, SpinStructure::RIGHT);
    BilinearTerm cc_right((Quark::CBAR), Quark::C, SpinStructure::RIGHT);
    
     FourQuarkOperator sduu = tensor_prod(sd_left, uu_left, ColorStructure::COLOR_DIAG);
    // Note we can use Fierz symmetry to tranform (sbar_a u_a)(ubar_b d_b) with
    // a color diagonal structure into (sbar_a d_b)(ubar_b u_a) with color
    // indices mixed
    FourQuarkOperator suud = tensor_prod(sd_left, uu_left, ColorStructure::COLOR_MIXED);
    FourQuarkOperator sddd = tensor_prod(sd_left, dd_left, ColorStructure::COLOR_DIAG);
    FourQuarkOperator sdss = tensor_prod(sd_left, ss_left, ColorStructure::COLOR_DIAG);
    FourQuarkOperator sdcc = tensor_prod(sd_left, cc_left, ColorStructure::COLOR_DIAG);
    // See note above on Fierz symmetry
    FourQuarkOperator sccd = tensor_prod(sd_left, cc_left, ColorStructure::COLOR_MIXED);
    
    // Covers P_1'' through P_5'', taken directly from the last line of eq. 67
    // in Masaaki's notes
    Ppp[0] = 5 * sduu + 5 * suud - 2 * sddd - 2 * sdss - 1 * sdcc - 1 * sccd;
    Ppp[0] *= 1.0 / (2 * std::sqrt(15));

    Ppp[1] = -1 * sduu - 1 * suud - 2 * sddd - 2 * sdss + 5 * sdcc + 5 * sccd;
    Ppp[1] *= 1.0 / (2 * std::sqrt(15));

    Ppp[2] = 1 * sduu - 1 * suud - 1 * sdcc + 1 * sccd;
    Ppp[2] *= 1.0 / 2;
    
    Ppp[3] = 1 * sduu - 1 * suud + 1 * sdcc - 1 * sccd;
    Ppp[3] *= 1.0 / 2;

    Ppp[4] = 1 * sduu + 1 * suud + 2 * sddd + 2 * sdss + 1 * sdcc + 1 * sccd;
    Ppp[4] *= 1.0 / (2 * std::sqrt(3));
    
    std::vector<BilinearTerm> up_like_bilinears = {{uu_right, cc_right}};
    for (BilinearTerm &qq : up_like_bilinears) {
        Ppp[5] += tensor_prod(sd_left, qq, ColorStructure::COLOR_DIAG);
        Ppp[6] += tensor_prod(sd_left, qq, ColorStructure::COLOR_MIXED);
        Ppp[7] += tensor_prod(sd_left, qq, ColorStructure::COLOR_DIAG);
        Ppp[8] += tensor_prod(sd_left, qq, ColorStructure::COLOR_MIXED);
    }
    std::vector<BilinearTerm> down_like_bilinears = {{dd_right, ss_right}};
    for (BilinearTerm &qq : down_like_bilinears) {
        Ppp[5] += tensor_prod(sd_left, qq, ColorStructure::COLOR_DIAG);
        Ppp[6] += tensor_prod(sd_left, qq, ColorStructure::COLOR_MIXED);
        Ppp[7] -= tensor_prod(sd_left, qq, ColorStructure::COLOR_DIAG);
        Ppp[8] -= tensor_prod(sd_left, qq, ColorStructure::COLOR_MIXED);
    }
    for (int i = 5; i < 7; i++) {
        Ppp[i] *= 1.0 / 2;
    }

    return Ppp;
}

static std::vector<FourQuarkOperator> get_operators_masaaki() {
    std::vector<FourQuarkOperator> Ppp; // P''_i from Masaaki's notes
    Ppp.resize(9);

    BilinearTerm sd_left((Quark::SBAR), Quark::D, SpinStructure::LEFT);
    BilinearTerm uu_left((Quark::UBAR), Quark::U, SpinStructure::LEFT);
    BilinearTerm dd_left((Quark::DBAR), Quark::D, SpinStructure::LEFT);
    BilinearTerm ss_left((Quark::SBAR), Quark::S, SpinStructure::LEFT);
    BilinearTerm cc_left((Quark::CBAR), Quark::C, SpinStructure::LEFT);
    BilinearTerm uu_right((Quark::UBAR), Quark::U, SpinStructure::RIGHT);
    BilinearTerm dd_right((Quark::DBAR), Quark::D, SpinStructure::RIGHT);
    BilinearTerm ss_right((Quark::SBAR), Quark::S, SpinStructure::RIGHT);
    BilinearTerm cc_right((Quark::CBAR), Quark::C, SpinStructure::RIGHT);

    //// Left-Left operator basis

    // Conventions are from Masaaki's notes, eq. 67
    // https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/tomii/NPR/NPR.pdf
    FourQuarkOperator sduu = tensor_prod(sd_left, uu_left, ColorStructure::COLOR_DIAG);
    // Note we can use Fierz symmetry to tranform (sbar_a u_a)(ubar_b d_b) with
    // a color diagonal structure into (sbar_a d_b)(ubar_b u_a) with color
    // indices mixed
    FourQuarkOperator suud = tensor_prod(sd_left, uu_left, ColorStructure::COLOR_MIXED);
    FourQuarkOperator sddd = tensor_prod(sd_left, dd_left, ColorStructure::COLOR_DIAG);
    FourQuarkOperator sdss = tensor_prod(sd_left, ss_left, ColorStructure::COLOR_DIAG);
    FourQuarkOperator sdcc = tensor_prod(sd_left, cc_left, ColorStructure::COLOR_DIAG);
    // See note above on Fierz symmetry
    FourQuarkOperator sccd = tensor_prod(sd_left, cc_left, ColorStructure::COLOR_MIXED);

    // Covers P_1'' through P_5'', taken directly from the last line of eq. 67
    // in Masaaki's notes
    Ppp[0] = 5 * sduu + 5 * suud - 2 * sddd - 2 * sdss - 1 * sdcc - 1 * sccd;
    Ppp[0] *= 1.0 / (2 * std::sqrt(15));

    Ppp[1] = -1 * sduu - 1 * suud - 2 * sddd - 2 * sdss + 5 * sdcc + 5 * sccd;
    Ppp[1] *= 1.0 / (2 * std::sqrt(15));

    Ppp[2] = 1 * sduu - 1 * suud - 1 * sdcc + 1 * sccd;
    Ppp[2] *= 1.0 / 2;

    Ppp[3] = 1 * sduu - 1 * suud + 1 * sdcc - 1 * sccd;
    Ppp[3] *= 1.0 / 2;

    Ppp[4] = 1 * sduu + 1 * suud + 2 * sddd + 2 * sdss + 1 * sdcc + 1 * sccd;
    Ppp[4] *= 1.0 / (2 * std::sqrt(3));


    //// Left-Right operators
    std::vector<BilinearTerm> up_like_bilinears = {{uu_right, cc_right}};
    for (BilinearTerm &qq : up_like_bilinears) {
        Ppp[5] += tensor_prod(sd_left, qq, ColorStructure::COLOR_DIAG);
        Ppp[6] += tensor_prod(sd_left, qq, ColorStructure::COLOR_MIXED);
        Ppp[7] += tensor_prod(sd_left, qq, ColorStructure::COLOR_DIAG);
        Ppp[8] += tensor_prod(sd_left, qq, ColorStructure::COLOR_MIXED);
    }
    std::vector<BilinearTerm> down_like_bilinears = {{dd_right, ss_right}};
    for (BilinearTerm &qq : down_like_bilinears) {
        Ppp[5] += tensor_prod(sd_left, qq, ColorStructure::COLOR_DIAG);
        Ppp[6] += tensor_prod(sd_left, qq, ColorStructure::COLOR_MIXED);
        Ppp[7] -= tensor_prod(sd_left, qq, ColorStructure::COLOR_DIAG);
        Ppp[8] -= tensor_prod(sd_left, qq, ColorStructure::COLOR_MIXED);
    }
    for (int i = 5; i < 9; i++) {
        Ppp[i] *= 1.0 / 2;
    }

    return Ppp;
}

std::vector<FourQuarkOperator> get_operators(OperatorBasis basis) {
    switch (basis) {
        case OperatorBasis::GREG: return get_operators_greg();
        case OperatorBasis::MASAAKI: return get_operators_masaaki();
        //Change here 
        //Using this because masaaki_qslash only supports 3 operators right now
        case OperatorBasis::MASAAKI_QSLASH: return get_operators_masaaki_qslash();
        //Change ends here
    }
}

FourQuarkTerm make_external_state(Quark q1, Quark q2, Quark q3, Quark q4) {
    BilinearTerm b1(q1, q2, SpinStructure::UNCONTRACTED);
    BilinearTerm b2(q3, q4, SpinStructure::UNCONTRACTED);
    return tensor_prod(b1, b2, ColorStructure::UNCONTRACTED);
}


//How do external states work? I think it is enforcing flavour part of Christoph paper
std::vector<FourQuarkOperator> get_external_states(ProjectorBasis basis) {
    FourQuarkOperator ext1 = make_external_state(Quark::S, Quark::DBAR, Quark::U, Quark::UBAR);
    FourQuarkOperator ext1_charmed = make_external_state(Quark::S, Quark::DBAR, Quark::C, Quark::CBAR);
    FourQuarkOperator ext1_stranged = make_external_state(Quark::S, Quark ::DBAR, Quark::S, Quark::SBAR);


    if (basis == ProjectorBasis::GREG_QSLASH) {
        std::vector<FourQuarkOperator> ret;
        for (int i = 0; i < 7; i++)
            ret.push_back(ext1);
        return ret;
    }
    //Change here 
    if(basis == ProjectorBasis::MASAAKI_QSLASH){
    	std::vector<FourQuarkOperator> ret;
    	for(int i = 0; i < 9; i++){
    		if(i == 1){
    			ret.push_back(ext1_charmed);
    			continue;
    		}
    		if(i == 4){
    			ret.push_back(ext1_stranged);
    			continue;
    		}
    		ret.push_back(ext1);
    	}
    	return ret;
    }
	//Cheange ends here
    std::vector<Quark> quarks = {Quark::U, Quark::D, Quark::S};
        if (basis == ProjectorBasis::MASAAKI) {
        quarks.push_back(Quark::C);
    }
    FourQuarkOperator ext2;
    for (Quark &q : quarks) {
        ext2 += make_external_state(Quark::S, Quark::DBAR, q, antiparticle(q));
    }

    std::vector<FourQuarkOperator> ret;
    ret.resize(7);
    // E_1 = E_2 = E_4 = E_5 = ext1 (Greg's thesis eq. 7.26)
    ret[0] = ext1;
    ret[1] = ext1;
    ret[3] = ext1;
    ret[4] = ext1;
    // E_3 = E_6 = E_7 = ext2 (Greg's thesis eq. 7.26)
    ret[2] = ext2;
    ret[5] = ext2;
    ret[6] = ext2;

    if (basis == ProjectorBasis::MASAAKI) {
        ret.push_back(ext1_charmed);
        ret.push_back(ext1_charmed);
    }

    return ret;
}

std::vector<TwoQuarkOp> get_subtraction_operators() {
    std::vector<TwoQuarkOp> ret = {
        TwoQuarkOp::SCALAR, TwoQuarkOp::DSLASH_RIGHT, TwoQuarkOp::DSLASH_LEFT
    };
    return ret;
}
