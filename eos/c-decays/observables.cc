/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2023-2025 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/observable-impl.hh>
#include <eos/c-decays/dq-to-l-nu.hh>
#include <eos/c-decays/dstarq-to-l-nu.hh>
#include <eos/c-decays/d-to-psd-l-nu.hh>
#include <eos/c-decays/lambdac-to-lambda-l-nu.hh>
#include <eos/c-decays/lambdac-to-proton-l-l.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Leptonic D decays
    // {{{
    ObservableGroup
    make_dq_to_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $D_q^{(*)+}\to \ell^+\nu$ decays)",
            R"(The option "l" selects the charged lepton flavor.)",
            {
                make_observable("D->lnu::BR", R"(\mathcal{B}(D^+ \to \ell^+\nu))",
                        Unit::None(),
                        &DqToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q"_ok, "d" } }),

                make_observable("D^*->lnu::BR", R"(\mathcal{B}(D^{*+} \to \ell^+\nu))",
                        Unit::None(),
                        &DstarqToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q"_ok, "d" } }),

                make_observable("D_s->lnu::BR", R"(\mathcal{B}(D_s^+ \to \ell^+\nu))",
                        Unit::None(),
                        &DqToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),

                make_observable("D_s^*->lnu::BR", R"(\mathcal{B}(D_s^{*+} \to \ell^+\nu))",
                        Unit::None(),
                        &DstarqToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Semileptonic D -> P(seudoscalar) decays
    // {{{

    // D -> K l nu
    // {{{
    ObservableGroup
    make_d_to_k_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $D\to K \ell^+ \nu$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("D->Klnu::dBR/dq2", R"(d\mathcal{B}(D\to K\ell^+ \nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &DToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),

                make_observable("D->Klnu::BR", R"(\mathcal{B}(D\to K\ell^+ \nu))",
                        Unit::None(),
                        &DToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),

                make_observable("D->Klnu::width", R"(\Gamma(D\to K\ell^+ \nu))",
                        Unit::None(),
                        &DToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),

                make_observable("D->Klnu::P(q2_min,q2_max)", R"(P(D\to K\ell^+ \nu))",
                        Unit::None(),
                        &DToPseudoscalarLeptonNeutrino::integrated_pdf_q2,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),

                make_observable("D->Klnu::P(q2)", R"(dP(D\to K\ell^+ \nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &DToPseudoscalarLeptonNeutrino::differential_pdf_q2,
                        std::make_tuple("q2"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Lambda_c decays
    // {{{
    ObservableGroup
    make_lambdac_to_lambda_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_c \to \Lambda \ell^+ \nu$ decays)",
            R"(The option "l" selects the charged lepton flavor.)",
            {
                make_observable("Lambda_c->Lambdalnu::BR", R"(\mathcal{B}(\Lambda_c^+ \to \Lambda \ell^+ \nu))",
                        Unit::None(),
                        &LambdaCToLambdaLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{}),

                make_observable("Lambda_c->Lambdalnu::dBR/dq2", R"(d\mathcal{B}/dq^2(\Lambda_c^+ \to \Lambda \ell^+ \nu))",
                        Unit::InverseGeV2(),
                        &LambdaCToLambdaLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{}),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Lambda_c -> p l l decays
    // {{{
    ObservableGroup
    make_lambdac_to_proton_l_l_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_c \to p \ell^+ \ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavor.)",
            {
                make_observable("Lambda_c->Protonll::Gamma", R"(\Gamma(\Lambda_c^+ \to p \ell^+ \ell^-))",
                        Unit::GeV(),
                        &LambdaCToProtonLeptonLepton::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{}),
 
                make_expression_observable("Lambda_c->Protonll::BR", R"(\mathcal{B}(\Lambda_c^+ \to p \ell^+ \ell^-))",
                        Unit::None(),
                        R"(
                        <<Lambda_c->Protonll::Gamma>> * [[life_time::Lambda_c]] / [[QM::hbar]]
                        )"),

                make_observable("Lambda_c->Protonll::dBR/dq2", R"(d\mathcal{B}/dq^2(\Lambda_c^+ \to p \ell^+ \ell^-))",
                        Unit::InverseGeV2(),
                        &LambdaCToProtonLeptonLepton::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{}),

                make_expression_observable("Lambda_c->Protonll::BR1/BR2", R"(\frac{\mathcal{B}_{[q2min1,q2max1]}}{\mathcal{B}_{[q2min2,q2max2]}}(\Lambda_c\to p \mu^+\mu^-))",
                        Unit::None(),
                        R"(
                        <<Lambda_c->Protonll::BR>>[q2_min=>q2_min_1,q2_max=>q2_max_1] 
                        / 
                        <<Lambda_c->Protonll::BR>>[q2_min=>q2_min_2,q2_max=>q2_max_2]
                        )"),

                make_expression_observable("Lambda_c->Protonll::BR1+BR2", R"(\mathcal{B}_{[q2min1,q2max1]}(\Lambda_c\to p \mu^+\mu^-) + \mathcal{B}_{[q2min2,q2max2]} (\Lambda_c\to p \mu^+\mu^-))",
                        Unit::None(),
                        R"(
                        <<Lambda_c->Protonll::BR>>[q2_min=>q2_min_1,q2_max=>q2_max_1] 
                        + 
                        <<Lambda_c->Protonll::BR>>[q2_min=>q2_min_2,q2_max=>q2_max_2]
                        )"),

                make_expression_observable("Lambda_c->Protonll::(BR1+BR2)/BR3", R"(\frac{\mathcal{B}_{[q2min1,q2max1]}+\mathcal{B}_{[q2min2,q2max2]}}{\mathcal{B}_{[q2min3,q2max3]}}(\Lambda_c\to p \mu^+\mu^-))",
                        Unit::None(),
                        R"(
                        (<<Lambda_c->Protonll::BR>>[q2_min=>q2_min_1,q2_max=>q2_max_1] 
                        + 
                        <<Lambda_c->Protonll::BR>>[q2_min=>q2_min_2,q2_max=>q2_max_2] ) 
                        / 
                        <<Lambda_c->Protonll::BR>>[q2_min=>q2_min_3,q2_max=>q2_max_3]
                        )"),

                make_observable("Lambda_c->Protonll::FL(q2)", R"(F_{\mathrm{L}}^{\Lambda_c^+ \to p \ell^+ \ell^-}(q^2))",
                        Unit::None(),
                        &LambdaCToProtonLeptonLepton::differential_fzero,
                        std::make_tuple("q2"),
                        Options{}),
    
                make_observable("Lambda_c->Protonll::AFB(q2)", R"(A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}(q^2))",
                        Unit::None(),
                        &LambdaCToProtonLeptonLepton::differential_a_fb_leptonic,
                        std::make_tuple("q2"),
                        Options{}),

                make_observable("Lambda_c->Protonll::FL_num", R"(\Gamma \cdot \langle F_{\mathrm{L}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                        Unit::GeV(),
                        &LambdaCToProtonLeptonLepton::integrated_fL_num,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{}),

                make_observable("Lambda_c->Protonll::AFB_num", R"(\Gamma \cdot \langle A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                        Unit::GeV(),
                        &LambdaCToProtonLeptonLepton::integrated_a_fb_leptonic_num,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{}),

                make_expression_observable("Lambda_c->Protonll::FL", R"(\langle F_{\mathrm{L}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                        Unit::None(),
                        R"(
                        <<Lambda_c->Protonll::FL_num>> / <<Lambda_c->Protonll::Gamma>>
                        )"),

                make_expression_observable("Lambda_c->Protonll::AFB", R"(\langle A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                        Unit::None(),
                        R"(
                        <<Lambda_c->Protonll::AFB_num>> / <<Lambda_c->Protonll::Gamma>>
                        )"),

                make_expression_observable("Lambda_c->Protonll::SigmaAFB", R"(\Sigma\langle A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                        Unit::None(),
                        R"(
                        0.5 * (<<Lambda_c->Protonll::AFB;cp-conjugate=false>> + <<Lambda_c->Protonll::AFB;cp-conjugate=true>>) 
                        )"),

                make_expression_observable("Lambda_c->Protonll::DeltaAFB", R"(\Delta\langle A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                        Unit::None(),
                        R"(
                        0.5 * (<<Lambda_c->Protonll::AFB;cp-conjugate=false>> - <<Lambda_c->Protonll::AFB;cp-conjugate=true>>) 
                        )"),
            }
        );
 
        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_c_decays_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in (semi)leptonic $c$-hadron decays",
            "",
            {
                // D_q^+ -> l^+ nu
                make_dq_to_l_nu_group(),

                // D -> K l^+ nu
                make_d_to_k_l_nu_group(),

                // Lc -> L l^+ nu
                make_lambdac_to_lambda_l_nu_group(),

                // Lc -> p l^+ l^-
                make_lambdac_to_proton_l_l_group()
            }
        );

        return ObservableSection(imp);
    }
}
