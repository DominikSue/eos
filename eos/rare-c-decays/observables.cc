/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2026 Dominik Suelmann
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
#include <eos/rare-c-decays/lambda-c-to-proton-l-l.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Lambda_c -> p l l decays
    // {{{
    ObservableGroup
    make_lambdac_to_proton_l_l_group()
    {
        auto imp = new Implementation<ObservableGroup>(
                R"(Observables in $\Lambda_c \to p \ell^+ \ell^-$ decays)",
                R"(The option "l" selects the charged lepton flavor.)",
                {
                    make_observable("Lambda_c->protonll::Gamma",
                                    R"(\Gamma(\Lambda_c^+ \to p \ell^+ \ell^-))",
                                    Unit::GeV(),
                                    &LambdaCToProtonLeptonLepton::integrated_decay_width,
                                    std::make_tuple("q2_min", "q2_max"),
                                    Options{}),

                    make_expression_observable("Lambda_c->protonll::BR",
                                               R"(\mathcal{B}(\Lambda_c^+ \to p \ell^+ \ell^-))",
                                               Unit::None(),
                                               R"(
                        <<Lambda_c->protonll::Gamma>> * [[life_time::Lambda_c]] / [[QM::hbar]]
                        )"),

                    make_observable("Lambda_c->protonll::dBR/dq2",
                                    R"(d\mathcal{B}/dq^2(\Lambda_c^+ \to p \ell^+ \ell^-))",
                                    Unit::InverseGeV2(),
                                    &LambdaCToProtonLeptonLepton::differential_branching_ratio,
                                    std::make_tuple("q2"),
                                    Options{}),

                    make_expression_observable("Lambda_c->protonll::BR1/BR2",
                                               R"(\frac{\mathcal{B}_{[q2min1,q2max1]}}{\mathcal{B}_{[q2min2,q2max2]}}(\Lambda_c\to p \mu^+\mu^-))",
                                               Unit::None(),
                                               R"(
                        <<Lambda_c->protonll::BR>>[q2_min=>q2_min_1,q2_max=>q2_max_1]
                        /
                        <<Lambda_c->protonll::BR>>[q2_min=>q2_min_2,q2_max=>q2_max_2]
                        )"),

                    make_expression_observable("Lambda_c->protonll::BR1+BR2",
                                               R"(\mathcal{B}_{[q2min1,q2max1]}(\Lambda_c\to p \mu^+\mu^-) + \mathcal{B}_{[q2min2,q2max2]} (\Lambda_c\to p \mu^+\mu^-))",
                                               Unit::None(),
                                               R"(
                        <<Lambda_c->protonll::BR>>[q2_min=>q2_min_1,q2_max=>q2_max_1]
                        +
                        <<Lambda_c->protonll::BR>>[q2_min=>q2_min_2,q2_max=>q2_max_2]
                        )"),

                    make_expression_observable("Lambda_c->protonll::(BR1+BR2)/BR3",
                                               R"(\frac{\mathcal{B}_{[q2min1,q2max1]}+\mathcal{B}_{[q2min2,q2max2]}}{\mathcal{B}_{[q2min3,q2max3]}}(\Lambda_c\to p \mu^+\mu^-))",
                                               Unit::None(),
                                               R"(
                        (<<Lambda_c->protonll::BR>>[q2_min=>q2_min_1,q2_max=>q2_max_1]
                        +
                        <<Lambda_c->protonll::BR>>[q2_min=>q2_min_2,q2_max=>q2_max_2])
                        /
                        <<Lambda_c->protonll::BR>>[q2_min=>q2_min_3,q2_max=>q2_max_3]
                        )"),

                    make_observable("Lambda_c->protonll::FL(q2)",
                                    R"(F_{\mathrm{L}}^{\Lambda_c^+ \to p \ell^+ \ell^-}(q^2))",
                                    Unit::None(),
                                    &LambdaCToProtonLeptonLepton::differential_f_l,
                                    std::make_tuple("q2"),
                                    Options{}),

                    make_observable("Lambda_c->protonll::AFB(q2)",
                                    R"(A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}(q^2))",
                                    Unit::None(),
                                    &LambdaCToProtonLeptonLepton::differential_a_fb_leptonic,
                                    std::make_tuple("q2"),
                                    Options{}),

                    make_observable("Lambda_c->protonll::FL_num",
                                    R"(\Gamma \cdot \langle F_{\mathrm{L}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                                    Unit::GeV(),
                                    &LambdaCToProtonLeptonLepton::integrated_f_l_num,
                                    std::make_tuple("q2_min", "q2_max"),
                                    Options{}),

                    make_observable("Lambda_c->protonll::AFB_num",
                                    R"(\Gamma \cdot \langle A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                                    Unit::GeV(),
                                    &LambdaCToProtonLeptonLepton::integrated_a_fb_leptonic_num,
                                    std::make_tuple("q2_min", "q2_max"),
                                    Options{}),

                    make_expression_observable("Lambda_c->protonll::FL",
                                               R"(\langle F_{\mathrm{L}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                                               Unit::None(),
                                               R"(
                        <<Lambda_c->protonll::FL_num>> / <<Lambda_c->protonll::Gamma>>
                        )"),

                    make_expression_observable("Lambda_c->protonll::AFB",
                                               R"(\langle A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                                               Unit::None(),
                                               R"(
                        <<Lambda_c->protonll::AFB_num>> / <<Lambda_c->protonll::Gamma>>
                        )"),

                    make_expression_observable("Lambda_c->protonll::SigmaAFB",
                                               R"(\Sigma\langle A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                                               Unit::None(),
                                               R"(
                        0.5 * (<<Lambda_c->protonll::AFB;cp-conjugate=false>> + <<Lambda_c->protonll::AFB;cp-conjugate=true>>)
                        )"),

                    make_expression_observable("Lambda_c->protonll::DeltaAFB",
                                               R"(\Delta\langle A_{\mathrm{FB}}^{\Lambda_c^+ \to p \ell^+ \ell^-}\rangle)",
                                               Unit::None(),
                                               R"(
                        0.5 * (<<Lambda_c->protonll::AFB;cp-conjugate=false>> - <<Lambda_c->protonll::AFB;cp-conjugate=true>>)
                        )"),
                });

        return ObservableGroup(imp);
    }

    // }}}

    ObservableSection
    make_rare_c_decays_section()
    {
        auto imp = new Implementation<ObservableSection>("Observables in rare $c$-hadron decays",
                                                         "",
                                                         { // Lc -> proton l^+ l^-
                                                           make_lambdac_to_proton_l_l_group() });

        return ObservableSection(imp);
    }
} // namespace eos
