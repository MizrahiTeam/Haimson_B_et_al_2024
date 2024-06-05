%% PB-related codes
% The code follows the figures order

% Task difficulty modulation of neuronal activity in the auditory cortex 
% Haimson B, Gilday OD, Lavi-Rudel A, Sagi H and Mizrahi A1*
% 
%% include in PB:



% Include: 
% save('AllVarsBlocks','AllData','n_mice','n_field','n_blocks_e','n_blocks_h','Seq_e','Seq_h','dprime_e','dprime_h','dprime_e_half','dprime_h_half',...
% 'Hit_rate_e','Hit_rate_h','FA_rate_e','FA_rate_h','relative_contrib_half_all','Res_GLM_half_all','Stat_half_all','R2_half_all',...
% 'R2_partial_half_all','Arr_GLM_half_all','Hit_e_tr_all','Hit_h_tr_all','CR_e_tr_all','CR_h_tr_all','FA_e_tr_all','FA_h_tr_all','Miss_e_tr_all',...
% 'Miss_h_tr_all','Hit_e_avg','Hit_h_avg','CR_e_avg','CR_h_avg','FA_e_avg','FA_h_avg','Miss_e_avg','Miss_h_avg','Res_all_list','PassiveData','NoGo_e_avg2','NoGo_h_avg2',...
% 'Go_dff_e_avg2','Go_dff_h_avg2','Go_e_avg2','Go_h_avg2','NonLR_list','AUC_e_avg','AUC_h_avg','AUC_e_blocks_avg','AUC_h_blocks_avg','Ac_e_all_tr','Ac_h_all_tr',...
% 'Ac_h_no2_tr','Ac_h_2_trx','mask','curB_half_all','Res_all')

% save('AllVarsLearning','CR_f_ng','Sig_f_ng','DATA','Go_idx','NeuronNum','NoGo_idx','AlignedD_N','AlignedD_NoGO2','FixedGo')

% save('AllVars_HaimsonB')

% ADD functions:
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1e and S1-1a.

dprime_h(dprime_h == -inf) = 0;
dprime_h_half(dprime_h_half == -inf) = 0;
for m = 1:n_mice
    for f = 1:n_field(m)
        figure; hold all % based on the whole blocks
        if m ~= 5
            x = interleave(dprime_h(m,f,1:n_blocks_h(m,f)),dprime_e(m,f,1:n_blocks_e(m,f)));
            for i = 1:numel(x)
                rectangle('Position',[i-.5 0 1 5],'FaceColor',[0 0 1/(mod(i,2)+1) .2],'EdgeColor','none')
            end
        else
            x = interleave(dprime_e(m,f,1:n_blocks_e(m,f)),dprime_h(m,f,1:n_blocks_h(m,f)));
            for i = 1:numel(x)
                rectangle('Position',[i-.5 0 1 5],'FaceColor',[0 0 1/(mod(i+1,2)+1) .2],'EdgeColor','none')
            end
        end
        plot(x)
        
        figure; hold all % based on the second half of the blocks
        if m ~= 5
            x = interleave(dprime_h_half(m,f,1:n_blocks_h(m,f)),dprime_e_half(m,f,1:n_blocks_e(m,f)));
            for i = 1:numel(x)
                rectangle('Position',[i-.5 0 1 5],'FaceColor',[0 0 1/(mod(i,2)+1) .2],'EdgeColor','none')
            end
        else
            x = interleave(dprime_e_half(m,f,1:n_blocks_e(m,f)),dprime_h_half(m,f,1:n_blocks_h(m,f)));
            for i = 1:numel(x)
                rectangle('Position',[i-.5 0 1 5],'FaceColor',[0 0 1/(mod(i+1,2)+1) .2],'EdgeColor','none')
            end
        end
        plot(x)
    end
end

% 1f

Dprime_e_behav = [];
Dprime_h_behav = [];

for m = 1:n_mice
    for f = 1:n_field(m)
        x = interleave(dprime_h(m,f,:),dprime_e(m,f,:));
        x2 = x(find(x,1,'first'):find(x,1,'last'));
        Dprime_e_behav = [Dprime_e_behav,x2(2:2:end)];
        Dprime_h_behav = [Dprime_h_behav,x2(1:2:end)];        
    end
end

figure; hold all
boxplot(Dprime_e_behav,'position',1)
boxplot(Dprime_h_behav,'position',2)
xticks(1:2)
ylim([-.5 4.5])
xlim([0.5 2.5])

[h,p] = ttest2(Dprime_e_behav,Dprime_h_behav)

% 1g

dpr1_e = NaN(n_mice,max(n_field),max(n_blocks_e,[],'all'));
dpr2_e = NaN(n_mice,max(n_field),max(n_blocks_e,[],'all'));            
dpr1_h = NaN(n_mice,max(n_field),max(n_blocks_h,[],'all'));
dpr2_h = NaN(n_mice,max(n_field),max(n_blocks_h,[],'all'));  
dpr_all_e = NaN(n_mice,max(n_field),max(n_blocks_e,[],'all'));            
dpr_all_h = NaN(n_mice,max(n_field),max(n_blocks_h,[],'all'));
for m = 1:n_mice
    for f = 1:n_field(m)
        for b = 1:n_blocks_e(m,f)
            x = Seq_e{m,f,b};
            H = sum(x(1:floor(end/2)) == 1)/sum(x(1:floor(end/2)) < 3);
            F = sum(x(1:floor(end/2)) == 3)/sum(x(1:floor(end/2)) > 2);
            if H == 1
                n = sum(x(1:floor(end/2)) < 3);
                H = (n-0.5)/n;
            end
            if F == 0
                n = sum(x(1:floor(end/2)) > 2);
                F = 0.5/n;
            end
            dpr1_e(m,f,b) = norminv(H) - norminv(F); 
            H = sum(x(ceil(end/2):end) == 1)/sum(x(ceil(end/2):end) < 3);
            F = sum(x(ceil(end/2):end) == 3)/sum(x(ceil(end/2):end) > 2);
            if H == 1
                n = sum(x(ceil(end/2):end) < 3);
                H = (n-0.5)/n;
            end
            if F == 0
                n = sum(x(ceil(end/2):end) > 2);
                F = 0.5/n;
            end
            dpr2_e(m,f,b) = norminv(H) - norminv(F);
            
            H = sum(x == 1)/sum(x < 3);
            F = sum(x == 3)/sum(x > 2);
            if H == 1
                n = sum(x(ceil(end/2):end) < 3);
                H = (n-0.5)/n;
            end
            if F == 0
                n = sum(x(ceil(end/2):end) > 2);
                F = 0.5/n;
            end
            dpr_all_e(m,f,b) = norminv(H) - norminv(F); 
        end
        for b = 1:n_blocks_h(m,f)
            x = Seq_h{m,f,b};
            H = sum(x(1:floor(end/2)) == 1)/sum(x(1:floor(end/2)) < 3);
            F = sum(x(1:floor(end/2)) == 3)/sum(x(1:floor(end/2)) > 2);
            if H == 1
                n = sum(x(1:floor(end/2)) < 3);
                H = (n-0.5)/n;
            end
            if F == 0
                n = sum(x(1:floor(end/2)) > 2);
                F = 0.5/n;
            end
            dpr1_h(m,f,b) = norminv(H) - norminv(F); 
            H = sum(x(ceil(end/2):end) == 1)/sum(x(ceil(end/2):end) < 3);
            F = sum(x(ceil(end/2):end) == 3)/sum(x(ceil(end/2):end) > 2);
            if H == 1
                n = sum(x(ceil(end/2):end) < 3);
                H = (n-0.5)/n;
            end
            if F == 0
                n = sum(x(ceil(end/2):end) > 2);
                F = 0.5/n;
            end
            dpr2_h(m,f,b) = norminv(H) - norminv(F);
            
            H = sum(x == 1)/sum(x < 3);
            F = sum(x == 3)/sum(x > 2);
            if H == 1
                n = sum(x(ceil(end/2):end) < 3);
                H = (n-0.5)/n;
            end
            if F == 0
                n = sum(x(ceil(end/2):end) > 2);
                F = 0.5/n;
            end
            dpr_all_h(m,f,b) = norminv(H) - norminv(F); 
        end
    end
end

x = dpr1_e(:);
y = dpr2_e(:);
x(x == -inf) = 0;
y(y == -inf) = 0;

x2 = dpr1_h(:);
y2 = dpr2_h(:);
x2(x2 == -inf) = 0;
y2(y2 == -inf) = 0;

figure; hold all
boxplot(x,'position',1)
boxplot(x2,'position',2)
boxplot(y,'position',3)
boxplot(y2,'position',4)
xticks(1:4)
xlim([.5 4.5])

[h,p] = ttest2(x,x2)
[h,p] = ttest2(y,y2)


% S1-1 b-c

hit_e = [];
hit_h = [];

for m = 1:n_mice
    for f = 1:n_field(m)
        for b = 1:n_blocks_e(m,f)
            hit_e = [hit_e;Hit_rate_e(m,f,b)];
        end
        for b = 1:n_blocks_h(m,f)
            hit_h = [hit_h;Hit_rate_h(m,f,b)];
        end
    end
end

hit_e(hit_e == 0) = [];
hit_h(hit_h == 0) = [];

figure; hold all
scatter(ones(numel(hit_e),1).*(1+(rand(numel(hit_e),1)-0.5)/5),hit_e,'filled')
scatter(ones(numel(hit_h),1).*(1+(rand(numel(hit_h),1)-0.5)/5)+1,hit_h,'filled')
line([0.8 1.2],[mean(hit_e) mean(hit_e)],'linewidth',2)
line([1.8 2.2],[mean(hit_h) mean(hit_h)],'linewidth',2)
ylim([0 1])
xlim([.5 2.5])

[h,p] = ranksum(hit_e,hit_h)


fa_e = [];
fa_h = [];

for m = 1:n_mice
    for f = 1:n_field(m)
        for b = 1:n_blocks_e(m,f)
            fa_e = [fa_e;FA_rate_e(m,f,b)];
        end
        for b = 1:n_blocks_h(m,f)
            fa_h = [fa_h;FA_rate_h(m,f,b)];
        end
    end
end

figure; hold all
scatter(ones(numel(fa_e),1).*(1+(rand(numel(fa_e),1)-0.5)/5),fa_e,'filled')
scatter(ones(numel(fa_h),1).*(1+(rand(numel(fa_h),1)-0.5)/5)+1,fa_h,'filled')
line([0.8 1.2],[mean(fa_e) mean(fa_e)],'linewidth',2)
line([1.8 2.2],[mean(fa_h) mean(fa_h)],'linewidth',2)
ylim([0 1])
xlim([.5 2.5])

[h,p] = ranksum(fa_e,fa_h)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2b

figure
h = heatmap(cat(1,relative_contrib_half_all{:,:}));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels = nan(size(h.XDisplayData));
grid off
colormap copper

% 2d

for m = 1:n_mice
    for f = 1:n_field(m)
        for n = 1:size(Res_GLM_half_all{m,f},1)
            p = Stat_half_all{m,f,n}.p;
            if p(end) < 0.05/numel(p)
                Res_diff{m,f}(n) = 1;
            else
                Res_diff{m,f}(n) = 0;
            end  
        end
    end
end

x = cat(1,curB_half_all{:});
figure
histogram(x(cat(2,Res_diff{:}) & cat(2,Res_all{:}),end),50,'Normalization','pdf')
hold on
xline(0,'k--')

% S2a

r = cat(2,Res_all{:});
x = [cat(2,R2_half_all{:})',cat(2,R2_partial_half_all{:,:,1})',cat(2,R2_partial_half_all{:,:,2})',cat(2,R2_partial_half_all{:,:,3})',cat(2,R2_partial_half_all{:,:,4})',cat(2,R2_partial_half_all{:,:,5})',...
    cat(2,R2_partial_half_all{:,:,6})',cat(2,R2_partial_half_all{:,:,7})'];
[~,I] = sort(mean(x(r,:)));
figure; hold all
y = mean(x(r,:));
scatter(1:8,y(I),'filled')
errorbar(1:8,y(I),std(x(r,:))./sqrt(size(x(r,:),1)),std(x(r,:))./sqrt(size(x(r,:),1)),'bp')
xticklabels(I)

[p,~,stats] = anova1([cat(2,R2_half_all{:})',cat(2,R2_partial_half_all{:,:,1})',cat(2,R2_partial_half_all{:,:,2})',cat(2,R2_partial_half_all{:,:,3})',cat(2,R2_partial_half_all{:,:,4})',cat(2,R2_partial_half_all{:,:,5})',...
    cat(2,R2_partial_half_all{:,:,6})',cat(2,R2_partial_half_all{:,:,7})']);
[c,m,h,gnames] = multcompare(stats);

% S2b

x = cat(1,relative_contrib_half_all{:,:});
r = cat(2,Res_all{:});

figure; hold all
bar(mean(x(r,:)).*100)
errorbar(1:7,mean(x(r,:)).*100,[],std(x(r,:).*100),'kp')

% S2c

for m = 1:n_mice
    for f = 1:n_field(m)
        for n = 1:size(Res_GLM_half_all{m,f},1)
            X = cat(1,Arr_GLM_half_all{m,f,:})';
            X(isnan(X)) = 0;
            k = 1;
            for w = -13:13
                [curB_half_all_win{m,f,k}(n,:),dev_half_all_win{m,f,k}(n),Stat_half_all_win{m,f,k,n}] = glmfit(X,circshift(Res_GLM_half_all{m,f}(n,:)',w),'normal','constant','off');
                cur_Ypred_half_all_win{m,f,k}(n,:) = [X]*curB_half_all_win{m,f,k}(n,:)';

                R2_half_all_win{m,f,k}(n) = corr(cur_Ypred_half_all_win{m,f,k}(n,:)',circshift(Res_GLM_half_all{m,f}(n,:)',w)).^2;
                R2_half_all_win_org{m,f,k}(n) = corr(cur_Ypred_half_all_win{m,f,k}(n,:)',circshift(Res_GLM_half_all{m,f}(n,:)',w)).^2;
                k = k + 1;
            end
        end
    end
end

clear R
for k = 1:numel(-13:13)
   R(:,k) = cat(2,R2_half_all_win{:,:,k})';
end

x = 1:numel(-13:13);
figure; hold all
plot(mean(R),'k')
patch('XData',[x flip(x)] , 'YData',[mean(R) + std(R) flip(mean(R)-std(R))],'FaceColor' , [0 0 0] , 'facealpha' , 0.2)
xlim([1 27])
ylim([0 1])
xticks(1:27)
xticklabels(strsplit(num2str(-13:13)))

[p,~,stats] = anova1(R);
[c,m,h,gnames] = multcompare(stats);

% S2d

for m = 1:n_mice
    for f = 1:n_field(m)
        for n = 1:size(Res_GLM_half_all{m,f},1)
            X = cat(1,Arr_GLM_half_all{m,f,:})';
            X(isnan(X)) = 0;
            k = 1;
            for w = -104:104
                cur_Ypred_half_all_win2{m,f,k}(n,:) = [X]*curB_half_all{m,f}(n,:)';
                R2_half_all_win2{m,f,k}(n) = corr(cur_Ypred_half_all_win2{m,f,k}(n,:)',circshift(Res_GLM_half_all{m,f}(n,:)',w)).^2;
                k = k + 1;
            end
        end
    end
end

clear R2
for k = 1:numel(-104:104)
   R2(:,k) = cat(2,R2_half_all_win2{:,:,k})';
end
    

x = 1:numel(-104:104);
figure; hold all
plot(mean(R2),'k')
patch('XData',[x flip(x)] , 'YData',[mean(R2) + std(R2) flip(mean(R2)-std(R2))],'FaceColor' , [0 0 0] , 'facealpha' , 0.2)
xlim([1 209])
ylim([-0.05 1])
xticks(1:26:209)
xticklabels({'-2\pi','-1.5\pi','-\pi','-0.5\pi','0','.5\pi','\pi','1.5\pi','2\pi'})

[p,~,stats] = anova1(R2);
[c,m,h,gnames] = multcompare(stats);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3 b,c

k = 1;
Sig_changed_e_h_CR = [];
for m = 1:n_mice
    for f = 1:n_field(m)
        X = CR_e_tr_all{m,f};
        Y = CR_h_tr_all{m,f};
        X2 = cat(2,X{:});
        Y2 = cat(2,Y{:});
        T = [];
        for n = 1:size(X2,1)
           if ranksum(X2(n,:),Y2(n,:),'tail','right') < 0.05
               Sig_changed_e_h_CR(k) = 1;
           elseif ranksum(X2(n,:),Y2(n,:),'tail','left') < 0.05
               Sig_changed_e_h_CR(k) = 2;
           else
               Sig_changed_e_h_CR(k) = 0;
           end
           T = [T;Sig_changed_e_h_CR(k)];
           k = k + 1;   
        end
        Sig_changed_e_h2{m,f} = T';
    end
end

% plot
X = cat(1,CR_e_avg{:});
Y = cat(1,CR_h_avg{:});

figure; hold all
gscatter(X(Res_all_list), Y(Res_all_list),Sig_changed_e_h_CR(Res_all_list),'krr','..',10,'filled')
plot(-0.1:0.01:.15,-0.1:0.01:.15,'k--')
xlim([-0.1 0.15])
ylim([-0.1 0.15])

figure
boxplot([X(Res_all_list),Y(Res_all_list)],'width',0.1)
ylim([-0.08 0.15])

[h,p] = ttest(X(Res_all_list),Y(Res_all_list))

% 3 f,g

X = cat(1,PassiveData.Dff_avg_e{:});
Y = cat(1,PassiveData.Dff_avg_h{:});

figure; hold all
gscatter(X(PassiveData.Res_all), Y(PassiveData.Res_all),PassiveData.Sig_changed_e_h_nogo(PassiveData.Res_all),'krr','..',10,'filled')
plot(-0.1:0.01:.15,-0.1:0.01:.15,'k--')
xlim([-0.1 0.15])
ylim([-0.1 0.15])

figure
boxplot([X(PassiveData.Res_all),Y(PassiveData.Res_all)],'width',0.1)
ylim([-0.08 0.15])

[h,p] = ttest(X(PassiveData.Res_all),Y(PassiveData.Res_all))



% S3 a-b


k = 1;
Sig_changed_e_h_NoGo2 = [];
for m = 1:n_mice
    for f = 1:n_field(m)
        X = CR_e_tr_all{m,f};
        Xf = FA_e_tr_all{m,f};
        Y = CR_h_tr_all{m,f};
        Yf = FA_h_tr_all{m,f};
        X2 = [cat(2,X{:}),cat(2,Xf{:})];
        Y2 = [cat(2,Y{:}),cat(2,Yf{:})];
        T = [];
        for n = 1:size(X2,1)
           if ranksum(X2(n,:),Y2(n,:),'tail','right') < 0.05
               Sig_changed_e_h_NoGo2(k) = 1;
           elseif ranksum(X2(n,:),Y2(n,:),'tail','left') < 0.05
               Sig_changed_e_h_NoGo2(k) = 2;
           else
               Sig_changed_e_h_NoGo2(k) = 0;
           end
           T = [T;Sig_changed_e_h_NoGo2(k)];
           k = k + 1;   
        end
        Sig_changed_e_hNoGo2{m,f} = T';
    end
end

X = cat(1,NoGo_e_avg2{:});
Y = cat(1,NoGo_h_avg2{:});
x = cellfun(@transpose,Sig_changed_e_hNoGo2,'UniformOutput',false);
sig = cat(1,x{:});

figure
gscatter(X(Res_all_list& NonLR_list), Y(Res_all_list& NonLR_list),sig(Res_all_list& NonLR_list),'krr','..',10,'filled')
hold on
plot(-0.1:0.01:.15,-0.1:0.01:.15,'k--')
xlim([-0.1 0.15])
ylim([-0.1 0.15])

figure
boxplot([X(Res_all_list& NonLR_list),Y(Res_all_list& NonLR_list)])
ylim([-0.06 .14])
[h,p] = ttest(X(Res_all_list& NonLR_list),Y(Res_all_list& NonLR_list))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4 b,c

k = 1;
Sig_changed_e_h_Hit = [];
for m = 1:n_mice
    for f = 1:n_field(m)
        X = Hit_e_tr_all{m,f};
        Y = Hit_h_tr_all{m,f};
        X2 = cat(2,X{:});
        Y2 = cat(2,Y{:});
        T = [];
        for n = 1:size(X2,1)
           if ranksum(X2(n,:),Y2(n,:),'tail','right') < 0.05
               Sig_changed_e_h_Hit(k) = 1;
           elseif ranksum(X2(n,:),Y2(n,:),'tail','left') < 0.05
               Sig_changed_e_h_Hit(k) = 2;
           else
               Sig_changed_e_h_Hit(k) = 0;
           end
           T = [T;Sig_changed_e_h_Hit(k)];
           k = k + 1;   
        end
        Sig_changed_e_h_Hit2{m,f} = T';
    end
end

X = cat(1,Hit_e_avg{:});
Y = cat(1,Hit_h_avg{:});

figure; hold all
gscatter(X(Res_all_list& NonLR_list), Y(Res_all_list& NonLR_list),Sig_changed_e_h_Hit(Res_all_list& NonLR_list),'krr','..',10,'filled')
plot(-0.155:0.01:.42,-0.155:0.01:.42,'k--')
xlim([-0.155 0.42])
ylim([-0.155 0.42])

figure
boxplot([X(Res_all_list& NonLR_list),Y(Res_all_list& NonLR_list)])
ylim([-.155 .42])
[h,p] = ttest(X(Res_all_list& NonLR_list),Y(Res_all_list& NonLR_list))

% 4 f,g

X = cat(1,Go_dff_e_avg2{:});
Y = cat(1,Go_dff_h_avg2{:});

figure; hold all
gscatter(X(PassiveData.Res_all), Y(PassiveData.Res_all),PassiveData.Sig_changed_e_h_go(PassiveData.Res_all),'krr','..',10,'filled')
plot(-0.13:0.01:.42,-0.13:0.01:.42,'k--')
xlim([-0.13 0.42])
ylim([-0.13 0.42])

figure
boxplot([X(PassiveData.Res_all),Y(PassiveData.Res_all)],'width',0.1)
ylim([-0.03 0.22])

[h,p] = ttest(X(PassiveData.Res_all),Y(PassiveData.Res_all))

%  S4 a,b

k = 1;
Sig_changed_e_h_Go2 = [];
for m = 1:n_mice
    for f = 1:n_field(m)
        X = Hit_e_tr_all{m,f};
        Xf = Miss_e_tr_all{m,f};
        Y = Hit_h_tr_all{m,f};
        Yf = Miss_h_tr_all{m,f};
        X2 = [cat(2,X{:}),cat(2,Xf{:})];
        Y2 = [cat(2,Y{:}),cat(2,Yf{:})];
        T = [];
        for n = 1:size(X2,1)
           if ranksum(X2(n,:),Y2(n,:),'tail','right') < 0.05
               Sig_changed_e_h_Go2(k) = 1;
           elseif ranksum(X2(n,:),Y2(n,:),'tail','left') < 0.05
               Sig_changed_e_h_Go2(k) = 2;
           else
               Sig_changed_e_h_Go2(k) = 0;
           end
           T = [T;Sig_changed_e_h_Go2(k)];
           k = k + 1;   
        end
        Sig_changed_e_hGo2{m,f} = T';
    end
end

X = cat(1,Go_e_avg2{:});
Y = cat(1,Go_h_avg2{:});
x = cellfun(@transpose,Sig_changed_e_hGo2,'UniformOutput',false);
sig = cat(1,x{:});

figure
gscatter(X(Res_all_list& NonLR_list), Y(Res_all_list& NonLR_list),sig(Res_all_list& NonLR_list),'krr','..',10,'filled')
hold on
plot(-0.1:0.01:.23,-0.1:0.01:.23,'k--')
xlim([-0.1 0.23])
ylim([-0.1 0.23])

figure
boxplot([X(Res_all_list& NonLR_list),Y(Res_all_list& NonLR_list)])
ylim([-0.06 .14])
[h,p] = ttest(X(Res_all_list& NonLR_list),Y(Res_all_list& NonLR_list))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5a
x = cat(2,AUC_e_avg{:})';
y = cat(2,AUC_h_avg{:})';
figure      
boxplot([x(Res_all_list & NonLR_list),y(Res_all_list & NonLR_list)],'width',0.1)
ylim([.5 1])
[h,p] = ttest(x(Res_all_list & NonLR_list),y(Res_all_list & NonLR_list))

% 5c
x = cat(2,AUC_e_avg{:})';
y = cat(2,AUC_h_avg{:})';
sig = cat(2,Sig_changed_e_h2{:});
sig_hit = cat(2,Sig_changed_e_h_Hit2{:});
r = cat(2,Res_all{:});
type1 = 0;
type2 = 2;
figure; hold all      
boxplot(x(r & NonLR_list & sig_hit ==type1 & sig == type1),'position',1,'width',0.3)
boxplot(x(r & NonLR_list & sig_hit ==type2),'position',2,'width',0.3)
boxplot(x(r & NonLR_list & sig == type2),'position',3,'width',0.3)
boxplot(y(r & NonLR_list & sig_hit ==type1 & sig == type1),'position',4,'width',0.3)
boxplot(y(r & NonLR_list & sig_hit ==type2),'position',5,'width',0.3)
boxplot(y(r & NonLR_list & sig == type2),'position',6,'width',0.3)
ylim([.5 1])
xlim([.5 6.5])
xticks(1:6)


[h,p] = ttest(x(r & NonLR_list & sig_hit ==type1 & sig == type1),y(r & NonLR_list & sig_hit ==type1 & sig == type1))
[h,p] = ttest(x(r & NonLR_list & sig_hit ==type2),y(r & NonLR_list & sig_hit ==type2))
[h,p] = ttest(x(r & NonLR_list & sig == type2),y(r & NonLR_list & sig == type2))

[h,p] = ttest2(x(r & NonLR_list & sig_hit ==type1 & sig == type1),x(r & NonLR_list & sig_hit ==type2))
[h,p] = ttest2(x(r & NonLR_list & sig_hit ==type1 & sig == type1),x(r & NonLR_list & sig ==type2))
[h,p] = ttest2(x(r & NonLR_list & sig_hit ==type2),x(r & NonLR_list & sig ==type2))

[h,p] = ttest2(y(r & NonLR_list & sig_hit ==type1 & sig == type1),y(r & NonLR_list & sig_hit ==type2))
[h,p] = ttest2(y(r & NonLR_list & sig_hit ==type1 & sig == type1),y(r & NonLR_list & sig ==type2))
[h,p] = ttest2(y(r & NonLR_list & sig_hit ==type2),y(r & NonLR_list & sig ==type2))

% 5d

x = cat(2,AUC_e_blocks_avg{:})';
y = cat(2,AUC_h_blocks_avg{:})';
figure
boxplot([x(PassiveData.Res_all),y(PassiveData.Res_all)],'width',0.1)
ylim([.49 1])

[h,p] = ttest(x(PassiveData.Res_all),y(PassiveData.Res_all))


% 5e

x1 = cellfun(@(x) nanmean(x,2) ,Ac_e_all_tr,'UniformOutput' ,false);
x2 = cellfun(@(x) nanmean(x,2) ,Ac_h_all_tr,'UniformOutput' ,false);
y1 = cellfun(@(x) nanmean(x,2) ,Ac_h_no2_tr,'UniformOutput' ,false);
y2 = cellfun(@(x) nanmean(x,2) ,Ac_h_2_trx,'UniformOutput' ,false);

figure; hold all
boxplot(cat(1,x1{mask}),'position',1)
boxplot(cat(1,x2{mask}),'position',2)
boxplot(cat(1,y2{mask}),'position',3)
boxplot(cat(1,y1{mask}),'position',4)
plot(1,mean(cat(1,x1{mask})),'dg')
plot(2,mean(cat(1,x2{mask})),'dg')
plot(3,mean(cat(1,y2{mask})),'dg')
plot(4,mean(cat(1,y1{mask})),'dg')
xticks(1:4)
ylim([0 1])
xlim([.5 4.5])

[h,p] = ttest2(cat(1,x1{mask}),cat(1,x2{mask}))
[h,p] = ttest2(cat(1,x2{mask}),cat(1,y1{mask}))
[h,p] = ttest2(cat(1,y2{mask}),cat(1,y1{mask}))
[h,p] = ttest2(cat(1,x2{mask}),cat(1,y2{mask}))

% S5a
x = cat(2,AUC_e_avg{:})';
y = cat(2,AUC_h_avg{:})';
sig = cat(2,Sig_changed_e_h2{:});
sig_hit = cat(2,Sig_changed_e_h_Hit2{:});
r = cat(2,Res_all{:});
type1 = 0;
type2 = 2;

x1 = x(r & NonLR_list & sig_hit ==type1 & sig == type1);
x2 = x(r & NonLR_list & sig_hit ==type2 );
x3 = x(r & NonLR_list & sig == type2);
y1 = y(r & NonLR_list & sig_hit ==type1 & sig == type1);
y2 = y(r & NonLR_list & sig_hit ==type2);
y3 = y(r & NonLR_list & sig == type2);

figure; hold all      
boxplot(y1-x1,'position',1,'width',0.3)
boxplot(y2-x2,'position',2,'width',0.3)
boxplot(y3-x3,'position',3,'width',0.3)
yline(0,'k--')
ylim([-0.25 0.32])
xlim([.5 3.5])
xticks(1:3)

[p,~,stats] = anova1([y1-x1;y2-x2;y3-x3],[ones(numel(x1),1);ones(numel(x2),1).*2;ones(numel(x3),1).*3]);
[c,m,h,gnames] = multcompare(stats);

% S5b

x = cat(2,AUC_e_avg{:})';
y = cat(2,AUC_h_avg{:})';
sig = cat(2,Sig_changed_e_h2{:});
sig_hit = cat(2,Sig_changed_e_h_Hit2{:});
r = cat(2,Res_all{:});
type1 = 0;
type2 = 1;
figure; hold all      
boxplot(x(r & NonLR_list & sig_hit ==type1 & sig == type1),'position',1,'width',0.3)
boxplot(x(r & NonLR_list & sig_hit ==type2),'position',2,'width',0.3)
boxplot(x(r & NonLR_list & sig == type2),'position',3,'width',0.3)
boxplot(y(r & NonLR_list & sig_hit ==type1 & sig == type1),'position',4,'width',0.3)
boxplot(y(r & NonLR_list & sig_hit ==type2),'position',5,'width',0.3)
boxplot(y(r & NonLR_list & sig == type2),'position',6,'width',0.3)
ylim([.5 1])
xlim([.5 6.5])
xticks(1:6)


[h,p] = ttest(x(r & NonLR_list & sig_hit ==type1 & sig == type1),y(r & NonLR_list & sig_hit ==type1 & sig == type1))
[h,p] = ttest(x(r & NonLR_list & sig_hit ==type2),y(r & NonLR_list & sig_hit ==type2))
[h,p] = ttest(x(r & NonLR_list & sig == type2),y(r & NonLR_list & sig == type2))

[h,p] = ttest2(x(r & NonLR_list & sig_hit ==type1 & sig == type1),x(r & NonLR_list & sig_hit ==type2))
[h,p] = ttest2(x(r & NonLR_list & sig_hit ==type1 & sig == type1),x(r & NonLR_list & sig ==type2))
[h,p] = ttest2(x(r & NonLR_list & sig_hit ==type2),x(r & NonLR_list & sig ==type2))

[h,p] = ttest2(y(r & NonLR_list & sig_hit ==type1 & sig == type1),y(r & NonLR_list & sig_hit ==type2))
[h,p] = ttest2(y(r & NonLR_list & sig_hit ==type1 & sig == type1),y(r & NonLR_list & sig ==type2))
[h,p] = ttest2(y(r & NonLR_list & sig_hit ==type2),y(r & NonLR_list & sig ==type2))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6b,c

figure; hold all
gscatter(CR_f_ng(4,:),CR_f_ng(5,:),Sig_f_ng,[0.6 0.6 0.6; 1 0 0 ; 1 0 0],'..',12,'filled')
plot(-0.05:0.01:.3,-0.05:0.01:.3,'k--')
xlim([-0.05 0.3])
ylim([-0.05 0.3])

figure
boxplot([CR_f_ng(4,:)',CR_f_ng(5,:)'],'width',0.1)

[h,p] = ttest(CR_f_ng(4,:)',CR_f_ng(5,:)')

% 6g

a = mean(mean(AlignedD_N(:,:,:,:,8:15),5),4);
b = squeeze(mean(a,2));
baseline_avg = mean(b(1:2,1:5));
b_std = squeeze(std(a,[],2))./sqrt(size(a,2));
baseline_avg_std = mean(b_std(1:2,1:5));
colors_list = [0 0 0;0.5 0.5 0.5; 0 1 0;1 0 1;1 0 0];
figure
x= 1:5;
for s = 1:3
    plot(x,b(s,1:5),'Color',colors_list(s,:),'LineWidth',2)
    patch('XData',[x flip(x)] , 'YData',[b(s,1:5) + b_std(s,1:5) flip(b(s,1:5)-b_std(s,1:5))],'FaceColor' , colors_list(s,:) , 'facealpha' , 0.2)
    hold on
end
ylim([0 0.12])
xticks(1:5)
xticklabels({'-0.5 oct','-0.25 oct','Go','0.25 oct','0.5 oct'})
ylabel('Average population response');
title('Average population response around Go frequency')

for i = 1:5
    [h(i),p(i)] = ttest(a(i,:,3),a(i,:,2));
end

% 6j

% [AlignedD_N,ind_N,LowFTail_N,HighFTail_N] = AlignData(NormalizeData(DATA),Go_idx(:,3),NeuronNum);
% [AlignedD_NoGO2,ind_N,LowFTail_N,HighFTail_N] = AlignData(NormalizeData(DATA),NoGo_idx(:,3),NeuronNum);

a = mean(mean(AlignedD_N(:,:,:,:,8:15),5),4);
a2 = mean(mean(AlignedD_NoGO2(:,:,:,:,8:15),5),4);

figure
boxplot([a(3,:,3)',a(4,:,3)',a2(3,:,1)',a2(4,:,1)'],'width',0.2)

[h,p] = ttest(a(3,:,3)',a(4,:,3)')
[h,p] = ttest(a2(3,:,1)',a2(4,:,1)')


% S6d

fr1 = [4,2,5,7,5,7];
Res_decay1 = [];
for s = 4:5
    k = 1;
    for m = 1:6
        NN = (sum(NeuronNum(1:(m-1)))+1):sum(NeuronNum(1:m));
        n = k:numel(NN)+k-1;
        Res_decay1(s,n) = mean(DATA(s,NN,fr1(m),:,8:15),4:5);
        k = k + numel(NN);
    end
end

figure
boxplot(Res_decay1(4:5,:)')
ylim([-0.1 0.6])

[h,p] = ttest(Res_decay1(4,:),Res_decay1(5,:))


% S6e
TransientData = [];
for m = 1:6
    NN = (sum(NeuronNum(1:(m-1)))+1):sum(NeuronNum(1:m));
    for s = 3:5
        if s > 4
            if FixedGo(m)
                TransientData(s,NN) = mean(squeeze(DATA(s,NN,NoGo_idx(m,4),:,8:15)),2:3);
            else
                TransientData(s,NN) = mean(squeeze(DATA(s,NN,Go_idx(m,4),:,8:15)),2:3);
            end
        else
            if FixedGo(m)
                TransientData(s,NN) = mean(squeeze(DATA(s,NN,NoGo_idx(m,s),:,8:15)),2:3);
            else
                TransientData(s,NN) = mean(squeeze(DATA(s,NN,Go_idx(m,s),:,8:15)),2:3);
            end
        end
    end
end

figure
boxplot(TransientData(3:5,:)')
ylim([-0.06 0.25])

[h,p] = ttest2(TransientData(3,:),TransientData(5,:))
[h,p] = ttest2(TransientData(3,:),TransientData(4,:))
[h,p] = ttest2(TransientData(4,:),TransientData(5,:))
