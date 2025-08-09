/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2025 Dominik Suelmann
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

#include <eos/c-decays/lambdac-to-proton-l-l.hh>
#include <eos/maths/complex.hh>
#include <eos/observable.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <test/test.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class LambdaCToProtonLeptonLeptonTest : public TestCase
{
public:
    LambdaCToProtonLeptonLeptonTest() : TestCase("lambdac_to_proton_l_l_test")
    {
    }

    virtual void
    run() const
    {
        {
            Parameters p = Parameters::Defaults();

            p["Lambda_c->Proton::a_0_time^V@DM2016"] = +0.83558637190456;
            p["Lambda_c->Proton::a_1_time^V@DM2016"] = -2.5698251471323;
            p["Lambda_c->Proton::a_2_time^V@DM2016"] = +9.8729125261567;
            p["Lambda_c->Proton::a_0_long^V@DM2016"] = +0.83253055920653;
            p["Lambda_c->Proton::a_1_long^V@DM2016"] = -2.3309542260103;
            p["Lambda_c->Proton::a_2_long^V@DM2016"] = +8.4088373584501;
            p["Lambda_c->Proton::a_0_perp^V@DM2016"] = +1.3629101111921;
            p["Lambda_c->Proton::a_1_perp^V@DM2016"] = -1.6996176062383;
            p["Lambda_c->Proton::a_2_perp^V@DM2016"] = +0.70894108564687;
            p["Lambda_c->Proton::a_0_time^A@DM2016"] = +0.72755668638733;
            p["Lambda_c->Proton::a_1_time^A@DM2016"] = -0.96674724145661;
            p["Lambda_c->Proton::a_2_time^A@DM2016"] = +0.82644748408508;
            p["Lambda_c->Proton::a_0_long^A@DM2016"] = +0.68653490341925;
            p["Lambda_c->Proton::a_1_long^A@DM2016"] = -0.90303820170714;
            p["Lambda_c->Proton::a_2_long^A@DM2016"] = +2.249274327598;
            p["Lambda_c->Proton::a_0_perp^A@DM2016"] = +0.68653490341925;
            p["Lambda_c->Proton::a_1_perp^A@DM2016"] = -0.68345507203321;
            p["Lambda_c->Proton::a_2_perp^A@DM2016"] = +0.69589476713259;
            p["Lambda_c->Proton::a_0_long^T@DM2016"] = +1.111320720066;
            p["Lambda_c->Proton::a_1_long^T@DM2016"] = -0.68870038548072;
            p["Lambda_c->Proton::a_2_long^T@DM2016"] = -2.8436650308711;
            p["Lambda_c->Proton::a_0_perp^T@DM2016"] = +0.63353942203683;
            p["Lambda_c->Proton::a_1_perp^T@DM2016"] = -1.0352768906776;
            p["Lambda_c->Proton::a_2_perp^T@DM2016"] = +1.4209332847548;
            p["Lambda_c->Proton::a_0_long^T5@DM2016"] = +0.62562392794376;
            p["Lambda_c->Proton::a_1_long^T5@DM2016"] = -1.1924360577673;
            p["Lambda_c->Proton::a_2_long^T5@DM2016"] = +3.7297035832191;
            p["Lambda_c->Proton::a_0_perp^T5@DM2016"] = +0.62562392794376;
            p["Lambda_c->Proton::a_1_perp^T5@DM2016"] = -1.3886937865559;
            p["Lambda_c->Proton::a_2_perp^T5@DM2016"] = +4.2242261969563;

            p["Lambda_c->Proton::res_a_rho"] = +0.54;
            p["Lambda_c->Proton::res_a_omega"] = +0.074;
            p["Lambda_c->Proton::res_a_phi"] = +0.106;
            p["Lambda_c->Proton::res_delta_rho"] = +0.00;
            p["Lambda_c->Proton::res_delta_omega_m_rho"] = +M_PI;
            p["Lambda_c->Proton::res_delta_phi_m_rho"] = +M_PI; //+0.00;
            p["gamma::rho^0"] = 0.1474;

            Options oo{
                {"model", "WET"},
                {"form-factors", "DM2016"},
                {"l", "mu"}};

            LambdaCToProtonLeptonLepton d(p, oo);

            const double eps = 3e-2;

            TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(0.75), 2.75291755e-7, eps);

            // the full phase-space region for muon
            TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(0.959, 1.122), 3.02e-7, eps);
        }
    }
} lambdac_to_proton_l_l_test;
