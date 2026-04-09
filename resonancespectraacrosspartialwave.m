datal0 = load("TmatrixvBFieldB10-20_L=0.mat");
datal2 = load("TmatrixvBFieldB10-20_L=2.mat");

figure; 
loglog(datal0.BFields_gauss_array, datal0.l0lp0);
hold on;
loglog(datal0.BFields_gauss_array, datal0.l0lp2);
loglog(datal2.BFields_gauss_array, datal2.l2lp2)
loglog(datal2.BFields_gauss_array, datal2.l2lp4)
xlabel("BField (Gauss")
ylabel("|T|^2")
legend({"L = 0->0", "L = 0->2", "L = 2->2", "L = 2->4"})