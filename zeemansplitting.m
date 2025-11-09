% This notebook calculates the Zeeman energies 

setup 

% calculate Zeeman energies 
gj = 1.167053; % from theory (cite) 
maxmj = 12; 
minmj = -12; 

dBField = 0.5;
BFieldlo = 1.0;
BFieldhi = 20.0;
numBF = floor((BFieldhi - BFieldlo)/dBField+0.1)+1; 
BFields_gauss_array = linspace(BFieldlo,BFieldhi,numBF);
base_thresholds = zeros(numBF); 
BFields = BFields_gauss_array/b0;
ZeemanEn = zeros(numBF, 7,7); 
ThresholdEns = zeros(numBF, length(HBmat)); 
mArr = -12:2:12;  
% thresholds 
targetQN = [-10,-10,4,-4]; % choose the channel to be the baseline 
targetQN = [-12,-12,0,0];
%energy = 1.6.d-3/t0 % collision energy
energy = 1.6e-8 


for mi = 1:length(mArr)
    for mj = 1:length(mArr)
        % define BField 
        for iBF = 1:numBF
            
            BField = BFields(iBF);
            thresholds = BField*diag(HBmat);  
            full_idx = find(all(Angular_QN_ULF == targetQN, 2));
            base_thresholds(iBF) = thresholds(full_idx); 
            m1 = mArr(mi); 
            m2 = mArr(mj); 
            
            Ez = muB*gj*BField*(m1+m2); 
            ZeemanEn(iBF, mi, mj) = Ez; 
        end     
    end 
end 

for iBF = 1:numBF
    ThresholdEns(iBF,:) = diag(HBmat)*BFields(iBF) ; 
end



% plot the threshold energies
figure;
hold on
legend_labels = cell(length(HBmat), 1);
for i = 1:length(HBmat)
    
%     legend_labels(i) = sprintf('Val1: %.2f, Val2: %.2f, Val3: %.2f, Val4: %.2f', ...
%         Angular_QN_ULF(i, 1), Angular_QN_ULF(i, 2), Angular_QN_ULF(i, 3), Angular_QN_ULF(i, 4));
    if Angular_QN_ULF(i,1) == -8 && Angular_QN_ULF(i,2) == -8 
        plot(BFields*b0, ThresholdEns(:,i)*t0-base_thresholds*t0, 'LineWidth', 3)
    else
        plot(BFields*b0, ThresholdEns(:,i)*t0-base_thresholds*t0, 'LineWidth', 1)
    end 
end
plot(BFields*b0, ones(length(BFields))*energy*t0, 'LineWidth', 2, linestyle = '--');
xlabel('B [Gauss]') 
ylabel('Threshold Energy [Kelvin]')

% make unique data - test this to see if it works (I want to see energy as
% a function of quantum number for each unique quantum number pair) 

m1 = Angular_QN_ULF(:, 1);
m2 = Angular_QN_ULF(:, 2);

[uniquePairs, idx] = unique([m1, m2], 'rows', 'stable');

% Take the last column of thresholdEns (energy at last B)
%E_threshold = ThresholdEns(:, end);
threshold_unique = ThresholdEns(15,idx)*t0;

idxnum = 1:length(uniquePairs); 

% Combine into one matrix
QN_Energies = [idxnum', uniquePairs, threshold_unique'];



% plotting
% figure; 
% hold on
% for i = 1:length(mArr)
%     for j = 1:length(mArr)
%         plot(BFields*b0, ZeemanEn(:,i,j)*t0);  
%         
%     end
% end 
% plot(BFields*b0, ones(length(BFields))*energy*t0, linestyle = '--');
% xlabel('B [Gauss]') 
% ylabel('Zeeman Energy [Kelvin]') 
% 
% 
% 
% % Now plot the difference for a given threshold, say m1 = -6, m2 = -6
% figure; 
% hold on
% for i = 1:length(mArr)
%     for j = 1:length(mArr)
%         %plot(BFields*b0, ZeemanEn(:,i,j)*t0 - base_thresholds*t0);  
%         plot(BFields*b0, ZeemanEn(:,i,j)*t0 - ZeemanEn(:,1,1)*t0);
%         
%     end
% end 
% plot(BFields*b0, ones(length(BFields))*energy*t0, 'LineWidth', 2, linestyle = '--');
% xlabel('B [Gauss]') 
% ylabel('Zeeman Energy [Kelvin]') 
