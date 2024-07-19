%% Plot the cumulative skew of OriC seqeunces after trimming, compressing and segregating

function [] = SkewInOriC (~)

%% Variables
MaxLen = 1000;                              % Maximum length of OriC sequence

DomainsFileName = {'MitoNCBI', 'DoriC12.1_bacteria', 'DoriC12.1_archaea', 'DoriC12.1_plasmid'};  
Path1= 'E:\Partha\data\DoriC\';         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% replace the file name and path if required

Domain_name = {'Mitochondria', 'Bacteria', 'Archaea', 'Plasmid'}; 

for OrgIdx = 1: length(DomainsFileName)
    Organism = DomainsFileName{OrgIdx};
    
    %% Import Data
    
    Data = readtable(strcat([Path1, Organism, '.csv']));
    
    %% Filter Data
    
    OriLen = zeros(size(Data,1),1);         % calculate length of origin sequences (may not be available in the data)
    for i = 1: size(Data,1)
        OriLen(i) = length(Data.oric_seq{i});
    end
    Data.OriLen = OriLen;                   % Calculate and store origin lengths
    if strcmp(Organism, 'MitoNCBI')         % special length cutoff for mito, as the oric seqeunce are of 10s of length
        MinLen = 16;
    else
        MinLen = 32;                        % minimum length of OriC sequence
    end
    
    Data = Data(Data.OriLen < MaxLen, :);   % remove long seqeunces
    Data = Data(Data.OriLen > MinLen, :);   % remove short sequences
    Seqs = Data.oric_seq;
    
    %% Process Data
    
    AllSeqsRY = alphaToNumericRY (Seqs);    % Convert the alphabatic sequence to numeric sequence with RY grouping
    [ApproxCoeffRY, DetailCoeffRY] = WaveletTransform (AllSeqsRY, MinLen);
    AvgDetailCoeffRY = sum(DetailCoeffRY)/size(DetailCoeffRY,1);
    AvgApproxCoeffRY = sum(ApproxCoeffRY)/size(ApproxCoeffRY,1);% Check V shape or A shape manually
    
    AllSeqsMK = alphaToNumericMK (Seqs);    % Convert the alphabatic sequence to numeric sequence with MK grouping
    [ApprCoeffMK, ~] = WaveletTransform (AllSeqsMK, MinLen);
    
    CorrRY= corr(DetailCoeffRY', AvgDetailCoeffRY');            % Calculate correlation
    
    %% Segregate Data
    
    if strcmp(Organism, 'MitoNCBI') || strcmp(Organism, 'DoriC12.1_bacteria') % for domains haveing V shape 
        LogicalArrayRY = CorrRY>=0; LogicalArrayMK = CorrRY<0;  % segregate indices; for V shape
    else
        LogicalArrayRY = CorrRY<0; LogicalArrayMK = CorrRY>=0;  % segregate indices; for A shape
    end
    
    BestAsRY = zeros(size(ApproxCoeffRY)); BestAsMK = zeros(size(ApprCoeffMK));
    BestAsRY(LogicalArrayRY,:) = ApproxCoeffRY(LogicalArrayRY,:);  BestAsMK(LogicalArrayMK,:) = ApprCoeffMK(LogicalArrayMK,:);
    BestAs = BestAsRY + BestAsMK;
    BestAvgAs = sum(BestAs)/size(BestAs,1); % Normalize
    
    %% Plot Data
    
    subplot(1,length(DomainsFileName),+OrgIdx)
    plot(BestAvgAs, 'LineWidth', 2)
    title(Domain_name{OrgIdx});
    if strcmp(Organism, 'MitoNCBI')         % set x-axis limit
        xlim([0,10])
    else
        xlim([0,20])
    end
    
    txt = (sprintf('%i RY + %i MK', sum(LogicalArrayRY), sum(LogicalArrayMK)));
    text(0, min(BestAvgAs), txt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    if OrgIdx == 2
        xlabel('Compressed sequence coordinate')
    end
    if OrgIdx == 1
        ylabel('Cumulative skew')
    end
end
end


%%
function [Appr_coeff_stored, Detail_coeff_stored] = WaveletTransform (S, MinLen)

Appr_coeff_stored = zeros(size(S,1),MinLen/2);    % Store Approx. Coefficients    
Detail_coeff_stored = zeros(size(S,1),MinLen/2);  % Store Detail Coefficients

for i = 1: size(S,1)
    Ns = S{i};
    
    %% Trim the sequence 
    
    El = 2^(floor(log2(size(Ns,2)))); Em = size(Ns,2);
    Tlen =(Em-El)*.5;
    if floor(Tlen)~=Tlen                    % trim the sequences
        Tlen = floor(Tlen);
        Ns = Ns((Tlen+1): (end-Tlen-1));    % trim for even
    else
        Ns = Ns(Tlen+1: end- Tlen);         % trim for odd
    end
    SeqCumsum = cumsum(Ns,2);                       % Calculate cumulative skew
    
    %% Compress the sequence
    
    Dn = log2(length(SeqCumsum))-log2(MinLen/2);    % Level of Wavelet transform
    if Dn ~= 0
        [C, L] = wavedec(SeqCumsum,Dn, 'haar');
        Appr_seq = C(1:L(1));                           % Approx. Coefficients
        Appr_seq = Appr_seq/(2^Dn);                     % Normalize the skew
        Detail_seq = C(L(1)+1:2*L(2));                  % Detailed Coefficients
    else
        Appr_seq = SeqCumsum; Detail_seq = zeros(size(SeqCumsum));
    end
    
    %% Store the seqeunce 
    
    Appr_coeff_stored(i,1:length(Appr_seq)) = Appr_seq;   % Store Approx. Coefficients
    Detail_coeff_stored(i,1:length(Appr_seq)) = Detail_seq;% Store Detail Coefficients
end

end



%% convert the set(cell) of alphabatic seq to numeric with purine-pyrimidine grouping

function [NSeqs] = alphaToNumericRY (Seqs)
NSeqs = cell(size(Seqs,1), 1);
for i = 1: size(Seqs,1)
    s = Seqs{i};
    s = char(lower(s));
    % Convert to numeric sequence
    Ns = zeros(size(s));
    Ns(s=='a') = 1;  Ns(s=='t') = -1;
    Ns(s=='g') = 1;  Ns(s=='c') = -1;
    NSeqs {i} = Ns;
end
end

function [NSeqs] = alphaToNumericMK (Seqs)
NSeqs = cell(size(Seqs,1), 1);
for i = 1: size(Seqs,1)
    s = Seqs{i};
    s = char(lower(s));
    % Convert to numeric sequence
    Ns = zeros(size(s));
    Ns(s=='a') = -1;  Ns(s=='t') = 1;
    Ns(s=='g') = 1;  Ns(s=='c') = -1;
    NSeqs {i} = Ns;
end
end









