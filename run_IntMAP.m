% Run IntMAP

clc
clear

% Processed PAS-Seq KO data file name or directory;
KO = 'PASseq_KO';
% Processed PAS-seq WT data file name or directory;
WT = 'PASseq_WT';
% Processed RNA-seq data file name or directory;
RNAseq = 'RNAseq';
% tuning parameter lambda
lambda = 10;
% p-value cutoff
pvalue_cutoff = 0.005;

load(WT);
Peak_WT = double(Peak);

Chr_WT = Chr;
clear C_flag Peak Chr gene

load(KO);
Peak_KO = double(Peak);

Chr_KO = Chr;
clear C_flag Peak Chr gene

Peak_Cutoff = 3;
n = 26569;
Loc_Cutoff = 40;


UniqueChr = intersect(unique(Chr_KO),unique(Chr_WT));

load mm10_uniqueTranscript CdsEnd CdsStart TxStart TxEnd ChrName Trend GeneName TranscriptName;

UniqueGene = unique(GeneName(1:n));

Candidate_Gene = cell(length(UniqueGene),2);
UTR_Start = zeros(length(UniqueGene),1);
UTR_End = zeros(length(UniqueGene),1);
UTR_Trend = zeros(length(UniqueGene),1);
UTR_Chr = cell(length(UniqueGene),1);
for i = 1:length(UniqueGene)
    ix = find(strcmp(GeneName(1:n),UniqueGene{i,1}));
    
    if length(unique(ChrName(ix)))==1
        if Trend{ix(1),1} == '+'
            
            [UTR_E, idx] = max(TxEnd(ix));
            UTR_S = min(CdsEnd(ix(idx)));
            UTR_E = UTR_E + 50;
            ix_WT = find(strcmp(Chr_WT,unique(ChrName(ix)))==1&Peak_WT(:,3)==1&Peak_WT(:,1)<=UTR_E&Peak_WT(:,1)>=UTR_S&Peak_WT(:,2)>=Peak_Cutoff);
            ix_KO = find(strcmp(Chr_KO,unique(ChrName(ix)))==1&Peak_KO(:,3)==1&Peak_KO(:,1)<=UTR_E&Peak_KO(:,1)>=UTR_S&Peak_KO(:,2)>=Peak_Cutoff);
            
            UTR_Start(i) = UTR_S;
            UTR_End(i) = UTR_E;
            UTR_Trend(i) = 1;
            UTR_Chr{i,1} = ChrName{ix(1),1};
            candidate = cat(1,Peak_WT(ix_WT,1),Peak_KO(ix_KO,1));
            
            if ~isempty(candidate)
                candidate = sort(unique(candidate));
            end
            
            if length(candidate)>1
                tmp = zeros(length(candidate),1);
                for j = 1:length(candidate)-1
                    if candidate(j+1) - candidate(j) <= Loc_Cutoff
                        if tmp(j) == 0 
                            tmp(j) = 0;
                            tmp(j+1) = 1;
                        else
                            tmp(j+1) = 0;
                        end
                    else
                        tmp(j) = 1;
                    end
                end
                if candidate(end) - candidate(end-1) <= Loc_Cutoff
                    if tmp(end-1) == 1 
                        tmp(end) = 0;
                    else
                        tmp(end) = 1;
                    end
                else
                    tmp(end) = 1;
                end
                candidate = candidate(tmp==1);
            end
            if length(candidate)>1
                tmp_KO = zeros(length(candidate),2);
                tmp_WT = zeros(length(candidate),2);
                tmp_KO(:,1) = candidate;
                tmp_WT(:,1) = candidate;
                
                for j = 1:length(candidate)
                    ix1 = find(abs(Peak_KO(ix_KO,1)-candidate(j))<=Loc_Cutoff);
                    ix2 = find(abs(Peak_WT(ix_WT,1)-candidate(j))<=Loc_Cutoff);
                    if ~isempty(ix1)
                        tmp_KO(j,2) = sum(Peak_KO(ix_KO(ix1),2));
                    end
                    if ~isempty(ix2)
                        tmp_WT(j,2) = sum(Peak_WT(ix_WT(ix2),2));
                    end
                end
                Candidate_Gene{i,1} = tmp_KO;
                Candidate_Gene{i,2} = tmp_WT;
            end
        else
    
            [UTR_S, idx] = min(TxStart(ix));
            UTR_E = max(CdsStart(ix(idx)));
            UTR_S = UTR_S - 50; 
            ix_WT = find(strcmp(Chr_WT,unique(ChrName(ix)))==1&Peak_WT(:,3)==0&Peak_WT(:,1)<=UTR_E&Peak_WT(:,1)>=UTR_S&Peak_WT(:,2)>=Peak_Cutoff);
            ix_KO = find(strcmp(Chr_KO,unique(ChrName(ix)))==1&Peak_KO(:,3)==0&Peak_KO(:,1)<=UTR_E&Peak_KO(:,1)>=UTR_S&Peak_KO(:,2)>=Peak_Cutoff);
            
            UTR_Start(i) = UTR_S;
            UTR_End(i) = UTR_E;
            UTR_Trend(i) = 0;
            UTR_Chr{i,1} = ChrName{ix(1),1};
            
            candidate = cat(1,Peak_WT(ix_WT,1),Peak_KO(ix_KO,1));
            
            if ~isempty(candidate)
                candidate = sort(unique(candidate));
            end
            
            if length(candidate)>1
                tmp = zeros(length(candidate),1);
                for j = 1:length(candidate)-1
                    if candidate(j+1) - candidate(j) <= Loc_Cutoff
                        if tmp(j) == 0 
                            tmp(j) = 0;
                            tmp(j+1) = 1;
                        else
                            tmp(j+1) = 0;
                        end
                    else
                        tmp(j) = 1;
                    end
                end
                if candidate(end) - candidate(end-1) <= Loc_Cutoff
                    if tmp(end-1) == 1 
                        tmp(end) = 0;
                    else
                        tmp(end) = 1;
                    end
                else
                    tmp(end) = 1;
                end
                candidate = candidate(tmp==1);
            end
            if length(candidate)>1
                tmp_KO = zeros(length(candidate),2);
                tmp_WT = zeros(length(candidate),2);
                tmp_KO(:,1) = candidate;
                tmp_WT(:,1) = candidate;
                
                for j = 1:length(candidate)
                    ix1 = find(abs(Peak_KO(ix_KO,1)-candidate(j))<=Loc_Cutoff);
                    ix2 = find(abs(Peak_WT(ix_WT,1)-candidate(j))<=Loc_Cutoff);
                    if ~isempty(ix1)
                        tmp_KO(j,2) = sum(Peak_KO(ix_KO(ix1),2));
                    end
                    if ~isempty(ix2)
                        tmp_WT(j,2) = sum(Peak_WT(ix_WT(ix2),2));
                    end
                end
                Candidate_Gene{i,1} = tmp_KO;
                Candidate_Gene{i,2} = tmp_WT;
            end
            
        end
    end
end

load(RNAseq);

P_KO = cell(length(UniqueGene),1); % relative proportion
Expression_KO = cell(length(UniqueGene),1); % expression level


P_WT = cell(length(UniqueGene),1);
Expression_WT = cell(length(UniqueGene),1);

maxiter = 100;
FailedIDX = zeros(length(UniqueGene),2);

for i = 1:length(UniqueGene)
    disp(['gene ',num2str(i),'is finished'])
    if mod(i,100)==0
        save IntMAP_Expression P_KO P_WT Expression_KO Expression_WT UniqueGene; % give the relative proportion and expression value for each transcripts in each gene in KO and WT
    end
    if ~isempty(DataMatrix_WT{i,1})
         try
         [P_KO{i,1}, Expression_KO{i,1}] = EM(DataMatrix_KO{i,1},PosInfo{i,1},maxiter,Candidate_Gene{i,1},lambda,UTR_Trend(i));
         catch
             FailedIDX(i,1)=1;
         end
        
         try
         [P_WT{i,1}, Expression_WT{i,1}] = EM(DataMatrix_WT{i,1},PosInfo{i,1},maxiter,Candidate_Gene{i,2},lambda,UTR_Trend(i));
         catch
             FailedIDX(i,2)=1;
         end
    end
end


pvalues = cell(length(P_KO),1);
R = zeros(length(P_KO),3);
pvalues_gene = -10*ones(length(P_KO),2);
for i = 1:length(P_KO)
    if ~isempty(Expression_KO{i,1})&&~isempty(Expression_WT{i,1})
        nn = length(Expression_KO{i,1});
        ptmp = -10*ones(nn-1,2);
        rtmp = zeros(nn-1,2);
        N1 = sum(Expression_KO{i,1})*10;
        N2 = sum(Expression_WT{i,1})*10;
        for j = 1:nn-1

                n1 = Expression_KO{i,1}(j)*10;
                n2 = Expression_WT{i,1}(j)*10;

            
            p0 = (n1+n2) / (N1+N2);
            
            % Expected counts under H0 (null hypothesis)
            n10 = N1 * p0;
            n20 = N2 * p0;

            % Chi-square test, by hand
            observed = [n1 N1-n1 n2 N2-n2];
            expected = [n10 N1-n10 n20 N2-n20];

            if n1/N1>= n2/N2        
                chi2stat = sum((observed-expected).^2 ./ expected);
                ptmp(j,1) = 1 - chi2cdf(chi2stat,1);
            else
                
                chi2stat = sum((observed-expected).^2 ./ expected);
                ptmp(j,2) = 1 - chi2cdf(chi2stat,1);
            end
            rtmp(j,1) = n1/N1;
            rtmp(j,2) = n2/N2;
            
        end
        pvalues{i,1} = ptmp;
        ix1 = find(ptmp(:,1)>=0&ptmp(:,1)<0.005);
        [a,im1] = min(ptmp(ix1,1));
        ix2 = find(ptmp(:,2)>=0&ptmp(:,2)<0.005);
        [b,im2]= min(ptmp(ix2,2));
        if ~isempty(a)&&~isempty(b)&&a<=b
            pvalues_gene(i,1) = a;
            R(i,1) = rtmp(ix1(im1),1);
            R(i,2) = rtmp(ix1(im1),2);
        elseif ~isempty(a)&&~isempty(b)&&a>b
            pvalues_gene(i,2) = b;
            R(i,1) = rtmp(ix2(im2),1);
            R(i,2) = rtmp(ix2(im2),2);
        elseif ~isempty(a)&&isempty(b)
            pvalues_gene(i,1) = a;
            R(i,1) = rtmp(ix1(im1),1);
            R(i,2) = rtmp(ix1(im1),2);
        elseif isempty(a)&&~isempty(b)
            pvalues_gene(i,2) = b;
            R(i,1) = rtmp(ix2(im2),1);
            R(i,2) = rtmp(ix2(im2),2);
        else
            R(i,1) = mean(rtmp(:,1));
            R(i,2) = mean(rtmp(:,2));
        end
        R(i,3) = 1;    
    end
end
idx1 = find(pvalues_gene(:,1)>=0&pvalues_gene(:,1)<pvalue_cutoff);
idx2 = find(pvalues_gene(:,2)>=0&pvalues_gene(:,2)<pvalue_cutoff);

G_KO = UniqueGene(idx1); % significant 3'UTR-APA gene in KO
G_WT = UniqueGene(idx2); % significant 3'UTR-APA gene in WT

Pvalue_KO = pvalues_gene(idx1,1); % pvalues for the significant genes in KO
Pvalue_WT = pvalues_gene(idx2,2); % pvalues for the significant genes in WT

save IntMAP_Results G_KO G_WT Pvalue_KO Pvalue_WT;
