% RM_04_FeatureViz.m
% 연도별 필드 dQ ratio 피처 시각화 + Lab 학습 범위 비교
clear; clc; close all;

%% Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir,'ExperimentalData','FeatureEngineering','Lab_RPT_Analysis','ver0317');
rawDir  = fullfile(projDir,'Rack_raw2mat');
rmDir   = fullfile(verDir,'RatioModel');
visDir  = fullfile(rmDir,'Visualization');
Q_nom=64; Np=2; num_segs=12;

%% Lab feature range
d=load(fullfile(verDir,'FeatureMatrix_ver0317.mat')); FM=d.FM;
mr=load(fullfile(verDir,'MasterRulers_PseudoOCV.mat'));
V_chg=mr.MasterRuler_ver0317.V_bounds_chg; V_dch=mr.MasterRuler_ver0317.V_bounds_dch;

chg_cols=arrayfun(@(i)sprintf('dQ_c_%02d',i),3:12,'UniformOutput',false);
dch_cols=arrayfun(@(i)sprintf('dQ_d_%02d',i),1:11,'UniformOutput',false);
dQ_c_r=FM{:,chg_cols}./sum(FM{:,chg_cols},2,'omitnan');
dQ_d_r=FM{:,dch_cols}./sum(FM{:,dch_cols},2,'omitnan');

lab_c_min=min(dQ_c_r,[],1,'omitnan'); lab_c_max=max(dQ_c_r,[],1,'omitnan');
lab_c_med=median(dQ_c_r,1,'omitnan');
lab_d_min=min(dQ_d_r,[],1,'omitnan'); lab_d_max=max(dQ_d_r,[],1,'omitnan');
lab_d_med=median(dQ_d_r,1,'omitnan');

%% Field data extraction
dataFiles={
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'),'Y2021','old',datetime(2021,6,3),[600 300];
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'),'Y2023','new',datetime(2023,10,16),[300 150];
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'),'Y2024','new',datetime(2024,9,9),[300 300];
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'),'Y2025','new',datetime(2025,7,11),[300 150];
};

years={'Y2021','Y2023','Y2024','Y2025'};
yr_colors={[0.22 0.43 0.74],[0.85 0.40 0.30],[0.29 0.62 0.35],[0.95 0.70 0.25]};

DC_r=nan(4,10); DD_r=nan(4,11); % dQ charge ratio (10 segs: 3~12), dQ dch (11 segs: 1~11)

for k=1:size(dataFiles,1)
    fpath=dataFiles{k,1}; dtype=dataFiles{k,3}; bd=dataFiles{k,4}; mdur=dataFiles{k,5};
    if ~exist(fpath,'file'), continue; end
    S=load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'),t=datetime(D.Time);
        elseif isduration(D.Date_Time),t=bd+D.Date_Time; else,t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
    else
        D=S.Raw;
        if isduration(D.Date_Time),t=bd+D.Date_Time; else,t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
    end
    tsec=seconds(t-t(1)); I_cell=I_rack/Np; thr=Q_nom*0.05/Np;
    chgS=local_segs(I_cell>thr); if ~isempty(chgS),dur=chgS(:,2)-chgS(:,1)+1;chgS=chgS(dur>=mdur(1),:);end
    dchS=local_segs(I_cell<-thr); if ~isempty(dchS),dur=dchS(:,2)-dchS(:,1)+1;dchS=dchS(dur>=mdur(2),:);end
    if isempty(chgS), continue; end

    [~,bi]=max(chgS(:,2)-chgS(:,1));
    v_c=V_avg(chgS(bi,1):chgS(bi,2)); i_c=I_cell(chgS(bi,1):chgS(bi,2)); t_c=tsec(chgS(bi,1):chgS(bi,2));
    % raw V (no smoothing)
    q_c=cumtrapz(t_c,abs(i_c))/3600;
    vm=v_c(1);qm=q_c(1);
    for ii=2:length(v_c),if v_c(ii)>vm(end),vm(end+1)=v_c(ii);qm(end+1)=q_c(ii);end;end
    dQ_chg=nan(1,num_segs);
    if length(vm)>1,QR=interp1(vm,qm,V_chg,'linear',NaN);dQ_chg=abs(diff(QR));end
    dc=dQ_chg(3:12); DC_r(k,:)=dc/sum(dc,'omitnan');

    if ~isempty(dchS)
        [~,bi_d]=max(dchS(:,2)-dchS(:,1));
        v_d=V_avg(dchS(bi_d,1):dchS(bi_d,2)); i_d=I_cell(dchS(bi_d,1):dchS(bi_d,2)); t_d=tsec(dchS(bi_d,1):dchS(bi_d,2));
        % raw V (no smoothing)
        q_d=cumtrapz(t_d,abs(i_d))/3600;
        vd=v_d(1);qd=q_d(1);
        for ii=2:length(v_d),if v_d(ii)<vd(end),vd(end+1)=v_d(ii);qd(end+1)=q_d(ii);end;end
        va=flip(vd);qa=flip(qd);
        dQ_dch=nan(1,num_segs);
        if length(va)>1,QR_d=interp1(va,qa,V_dch,'linear',NaN);dQ_dch=abs(diff(QR_d));end
        dd=dQ_dch(1:11); DD_r(k,:)=dd/sum(dd,'omitnan');
    end
end

%% ===== Visualization =====
seg_c = 3:12;   % charge seg indices
seg_d = 1:11;   % discharge seg indices

fig=figure('Position',[50 50 1200 540],'Color','w');
sgtitle('Field dQ Ratio Feature Profiles per Year vs. Lab Training Range',...
    'FontSize',13,'FontWeight','bold');

% --- Charge ---
ax1=subplot(1,2,1);
hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'GridColor',[0.88 0.88 0.88],'GridAlpha',1,'FontSize',10,'TickDir','out');

% Lab range shading
fill(ax1,[seg_c fliplr(seg_c)],[lab_c_min fliplr(lab_c_max)],...
    [0.75 0.85 0.95],'FaceAlpha',0.4,'EdgeColor','none','DisplayName','Lab range');
% Lab median
plot(ax1,seg_c,lab_c_med,'b--','LineWidth',1.5,'DisplayName','Lab median');

% Field per year
for k=1:4
    dc=DC_r(k,:);
    valid=~isnan(dc);
    if any(valid)
        plot(ax1,seg_c(valid),dc(valid),'-o','Color',yr_colors{k},'LineWidth',2,...
            'MarkerSize',7,'MarkerFaceColor',yr_colors{k},'DisplayName',years{k});
        % NaN markers
        if any(~valid)
            plot(ax1,seg_c(~valid),zeros(1,sum(~valid)),'x','Color',yr_colors{k},...
                'MarkerSize',10,'LineWidth',2,'HandleVisibility','off');
        end
    end
end

xlabel(ax1,'Charge Segment Index','FontSize',11);
ylabel(ax1,'dQ Ratio','FontSize',11);
title(ax1,'Charge dQ Ratio (Seg 03–12)','FontSize',11,'FontWeight','bold');
legend(ax1,'Location','northwest','FontSize',9,'Box','off');
xlim(ax1,[2.5 12.5]); set(ax1,'XTick',seg_c);

% --- Discharge ---
ax2=subplot(1,2,2);
hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'GridColor',[0.88 0.88 0.88],'GridAlpha',1,'FontSize',10,'TickDir','out');

fill(ax2,[seg_d fliplr(seg_d)],[lab_d_min fliplr(lab_d_max)],...
    [0.75 0.85 0.95],'FaceAlpha',0.4,'EdgeColor','none','DisplayName','Lab range');
plot(ax2,seg_d,lab_d_med,'b--','LineWidth',1.5,'DisplayName','Lab median');

for k=1:4
    dd=DD_r(k,:);
    valid=~isnan(dd);
    if any(valid)
        plot(ax2,seg_d(valid),dd(valid),'-o','Color',yr_colors{k},'LineWidth',2,...
            'MarkerSize',7,'MarkerFaceColor',yr_colors{k},'DisplayName',years{k});
        if any(~valid)
            plot(ax2,seg_d(~valid),zeros(1,sum(~valid)),'x','Color',yr_colors{k},...
                'MarkerSize',10,'LineWidth',2,'HandleVisibility','off');
        end
    end
end

xlabel(ax2,'Discharge Segment Index','FontSize',11);
ylabel(ax2,'dQ Ratio','FontSize',11);
title(ax2,'Discharge dQ Ratio (Seg 01–11)','FontSize',11,'FontWeight','bold');
legend(ax2,'Location','northeast','FontSize',9,'Box','off');
xlim(ax2,[0.5 11.5]); set(ax2,'XTick',seg_d);

% Footer note
annotation(fig,'textbox',[0.02 0.01 0.96 0.04],'String',...
    '× : NaN (segment not available in field data) | Shaded: Lab training range (min–max) | Dashed: Lab median',...
    'EdgeColor','none','FontSize',9,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');

saveas(fig, fullfile(visDir,'RM04_Field_FeatureProfiles_RawV.png'));
fprintf('Saved: RM04_Field_FeatureProfiles_RawV.png\n');

%% Helpers
function segs=local_segs(mask)
    d=[0;diff(mask(:))]; starts=find(d==1); ends=find(d==-1)-1;
    if mask(1),starts=[1;starts];end; if mask(end),ends=[ends;length(mask)];end
    if isempty(starts)||isempty(ends),segs=[];return;end; segs=[starts,ends];
end
