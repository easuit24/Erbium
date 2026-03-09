%cd("Tmatanalysis_fullvsdipolevsBorn/")
born02 = load("born02.mat");
born22 = load("born22.mat");
born24 = load("born24.mat");

fullcalc02_off = load("fullcalc02_offresonance.mat");
fullcalc22_off = load("fullcalc22_offresonance.mat");
fullcalc24_off = load("fullcalc24_offresonance.mat");
dipole02_off = load("dipoleonly02_offresonance.mat");
dipole22_off = load("dipoleonly22_offresonance.mat");
dipole24_off = load("dipoleonly24_offresonance.mat");

fullcalc02_on = load("fullcalc02_onresonance.mat");
fullcalc22_on = load("fullcalc22_onresonance.mat");
fullcalc24_on = load("fullcalc24_onresonance.mat");
dipole02_on = load("dipoleonly02_onresonance.mat");
dipole22_on = load("dipoleonly22_onresonance.mat");
dipole24_on = load("dipoleonly24_onresonance.mat");

% make the relevant plots for all of these - one for on resonance and one
% for off resonance

figure; 
loglog(born02.energies2, born02.temparr, ":", LineWidth = 2)
hold on
loglog(born22.energies2, born22.temparr, ":", LineWidth = 2)
loglog(born24.energies2, born24.temparr, ":", LineWidth = 2)
loglog(fullcalc02_off.energies, fullcalc02_off.temparr, "-", LineWidth = 2)
loglog(fullcalc22_off.energies, fullcalc22_off.temparr, "-", LineWidth = 2)
loglog(fullcalc24_off.energies, fullcalc24_off.temparr, "-", LineWidth = 2)
loglog(dipole02_off.energies, dipole02_off.temparr, "--", LineWidth = 2)
loglog(dipole22_off.energies, dipole22_off.temparr, "--", LineWidth = 2)
loglog(dipole24_off.energies, dipole24_off.temparr, "--", LineWidth = 2)

xlabel("Energy (atomic units)")
ylabel("|T|^2")
title("Off Resonance")


figure; 
loglog(born02.energies2, born02.temparr, ":", LineWidth = 2)
hold on
loglog(born22.energies2, born22.temparr, ":", LineWidth = 2)
loglog(born24.energies2, born24.temparr, ":", LineWidth = 2)
loglog(fullcalc02_on.energies, fullcalc02_on.temparr, "-", LineWidth = 2)
loglog(fullcalc22_on.energies, fullcalc22_on.temparr, "-", LineWidth = 2)
loglog(fullcalc24_on.energies, fullcalc24_on.temparr, "-", LineWidth = 2)
loglog(dipole02_on.energies, dipole02_on.temparr, "--", LineWidth = 2)
loglog(dipole22_on.energies, dipole22_on.temparr, "--", LineWidth = 2)
loglog(dipole24_on.energies, dipole24_on.temparr, "--", LineWidth = 2)
xlabel("Energy (atomic units)")
ylabel("|T|^2")
title("On Resonance")