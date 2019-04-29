%--------------------------------------------------------------------------
%                        ECHELLE DIAGRAM + MODE ID
%
% Directly adapted from TESS_chaplin_first_targets.m
%--------------------------------------------------------------------------
function oscillation = plato_mode_ID_echelle(choix,data,oscillation)
%%
vert      = [.11 .62 .47];
bleu      = [0    0.25 0.53];
vermillon = [0.89 0.29 0.20];

%... ID variables from input structures
dir_plot    = choix.dir_plot;
nu          = data.nu;
PSD         = data.PSD;
PSD_w       = data.PSD_w;
PIC         = oscillation.ID_target;
num_target  = oscillation.num_target;
Dnu         = oscillation.Dnu;
Dnu_err     = oscillation.Dnu_err;
Teff        = oscillation.Teff;
nu_max      = oscillation.nu_max;
nu_max_err  = oscillation.nu_max_err;
nu_peak     = oscillation.psd_H0.nu_peak;
N_n         = oscillation.N_n;

clear nu_X nu_Y diag_ech nu_0_bis
oscillation.R_Rsun     = 0;

if oscillation.detection == 1
    
     %... Scaling relation
    [M_Msun,sig_M,R_Rsun,sig_R,logg,sig_logg,ro_ro_sun,sig_ro] = scaling_law_RG(nu_max,nu_max_err,Dnu,max([Dnu_err 0.01]),Teff,100);
    
    %... Echelle diagram associated with mean large spacing
    nu_inf_ech = oscillation.psd_H0.nu_range(1);
    nu_sup_ech = oscillation.psd_H0.nu_range(2);
    nn   = find(nu > nu_inf_ech & nu < nu_sup_ech);
    norma = 1;
    if nu_max < 500
        N_liss = 1;
    else
        N_liss = 3;
    end
    PSD_input_ech = lissage3(PSD_w(nn),N_liss);
    interpol = 1;
    [nu_X,nu_Y,diag_ech,nu_0_bis] = diag_echelle(PSD_input_ech,nu(nn),Dnu,norma,interpol,0);
    
    %... Detected peaks folded on the echelle diagram
    offset_X   = 0;
    freq_ech_X = mod(nu_peak - nu_inf_ech,Dnu) - offset_X;
    freq_ech_Y = nu_peak - freq_ech_X;
    
    figure
    fontzi = 14;
    imagesc(nu_X,nu_Y,-diag_ech,[-5 0])
    hold
    plot(freq_ech_X,freq_ech_Y,'s','color','b','linewidth',2)
    colormap gray
    set(gca,'YDir','normal','fontsize',fontzi,'ytick',round(nu_Y))
    
    %......................................................................
    %                  Comparison with universal pattern
    %......................................................................
    %... Revision of echelle diagram and coordinates on H0 peaks on the
    %    echelle diagram, according to new Dnu_0
    [nu_X,nu_Y,diag_ech,nu_0_bis] = diag_echelle(PSD_input_ech,nu(nn),Dnu,norma,interpol,0);
    max_val_ech = [];
    
    %... Universal pattern
    if Teff < 5800
        fprintf('==================== \n')
        fprintf('RG Universal pattern \n')
        fprintf('==================== \n')
        [nu_up_l0,nu_up_l1,nu_up_l2] = diag_echelle_UP_RG_V2(nu,PSD,Dnu,nu_max,PIC,dir_plot,N_liss,N_n,max_val_ech,offset_X,nu_inf_ech);
        freq_ech_X = mod(nu_peak - nu_inf_ech,Dnu) - offset_X;
        freq_ech_Y = nu_peak - freq_ech_X;
        plot(freq_ech_X,freq_ech_Y,'s','color','b','linewidth',2)
        title(['PIC ' PIC ': numax = ' num2str(nu_max,'%6.0f') ', M=' num2str(M_Msun,'%6.2f') ', R=' num2str(R_Rsun,'%6.2f') ', T=' num2str(Teff,'%6.0f')],'interpreter','latex','fontsize',fontzi+2)
        nom_figue = [dir_plot 'KIC_' PIC '_echelle_UP.png'];
        set(gcf,'PaperPositionMode','auto')
        eval(['print -dpng -loose ' nom_figue])
    end
    
    %... Universal pattern MS stars
    if Teff > 5800
        %... RG/SG
        %Aeps = 0.60; Beps = 0.52; Ceps = 0.005;
        %eps  = Aeps + log10(Dnu)*Beps + Dnu*Ceps;
        fprintf('==================== \n')
        fprintf('MS Universal pattern \n')
        fprintf('==================== \n')
        
        %... MS
        eps = 4.62 - 5.55e-4*Teff - 0.25; % Note: - 0.25 from Warrick Ball
                
        [nu_up_l0,nu_up_l1,nu_up_l2,nu_inf] = diag_echelle_UP_MS(nu,PSD_w,Dnu,nu_max,PIC,dir_plot,N_liss,N_n,max_val_ech,offset_X,nu_inf_ech,eps);
        freq_ech_X = mod(nu_peak - nu_inf,Dnu) - offset_X;
        freq_ech_Y = nu_peak - freq_ech_X;
        plot(freq_ech_X,freq_ech_Y,'s','color','b','linewidth',2)
        title(['PIC ' PIC ': numax = ' num2str(nu_max,'%6.0f') ', M=' num2str(M_Msun,'%6.2f') ', R=' num2str(R_Rsun,'%6.2f') ', T=' num2str(Teff,'%6.0f')],'interpreter','latex','fontsize',fontzi+2)

    end
    
    fig_number = get(gcf,'Number');
    
    %......................................................................
    %                      REVISED universal pattern
    %......................................................................
    %... Check shifts
    [shift_l0,shift_l1] = plato_mode_ID_shift_wrt_UP(nu,PSD,Dnu,nu_up_l0,PIC,dir_plot);
    
    %... Figure followed
    figure(fig_number)
    
    %... Revised epsilon
    alpha       = 2*0.57/(nu_max/Dnu)^2; % D'apres email Benoit 30.1.2019
    n_l0        = round(nu_up_l0/Dnu);%floor(nu_up_l0/Dnu);
    eps_revised = eps + shift_l0/Dnu;
    
    %... l=0
    nu_up_l0_revised = (n_l0 + eps_revised + alpha/2*(n_l0 - nu_max/Dnu).^2)*Dnu;
    freq_ech_X = mod(nu_up_l0_revised - nu_inf,Dnu) - offset_X;
    freq_ech_Y = nu_up_l0 + shift_l0 - freq_ech_X;
    plot(freq_ech_X,freq_ech_Y,'--','color',bleu,'linewidth',2)

    d01_revised  = shift_l1/Dnu;
    nu_up_l1_revised = (n_l0 + 1/2 + eps_revised - d01_revised + alpha/2*(n_l0 - nu_max/Dnu).^2)*Dnu; 
    freq_ech_X = mod(nu_up_l1_revised - nu_inf,Dnu) - offset_X;
    freq_ech_Y = nu_up_l1_revised - freq_ech_X;
    plot(freq_ech_X,freq_ech_Y,'--','color',vert,'linewidth',2)
    
    %... For l=2 i prefer keeping the "theoretical" d02 with revised
    %    epsilon because it often fails
    d02              = (0.131 -0.033*log10(Dnu))+ 0.03;
    %      d02_revised      = shift_l2/Dnu;
    nu_up_l2_revised = (n_l0 + 2/2 + eps_revised - d02 + alpha/2*(n_l0 - nu_max/Dnu).^2)*Dnu;
    freq_ech_X = mod(nu_up_l2_revised - nu_inf,Dnu) - offset_X;
    freq_ech_Y = nu_up_l2_revised - freq_ech_X;
    plot(freq_ech_X,freq_ech_Y,'--','color',vermillon,'linewidth',2)
    
    %......................................................................
    %                      Mode ID l = (0,1) and l=2
    %......................................................................
    %... l = (0,1)
    [nu_l0,nu_l1] = plato_mode_ID_echelle_subroutine_l0_l1(nu_peak,nu_up_l0_revised,nu_up_l1_revised,Dnu);
    freq_ech_X_l0 = mod(nu_l0 - nu_inf,Dnu) - offset_X;
    freq_ech_Y_l0 = nu_l0 - freq_ech_X_l0;
    freq_ech_X_l1 = mod(nu_l1 - nu_inf,Dnu) - offset_X;
    freq_ech_Y_l1 = nu_l1 - freq_ech_X_l1;
    
    figure(fig_number)
    plot(freq_ech_X_l0,freq_ech_Y_l0,'d','color',bleu,'markersize',8,'linewidth',2)
    plot(freq_ech_X_l1,freq_ech_Y_l1,'d','color',vert,'markersize',8,'linewidth',2)

    %... l = 2
    nu_l2         = plato_mode_ID_echelle_subroutine_l2(nu,PSD,nu_l0,nu_l1,nu_up_l2_revised,d02,oscillation);
    freq_ech_X_l2 = mod(nu_l2 - nu_inf,Dnu) - offset_X;
    freq_ech_Y_l2 = nu_l2 - freq_ech_X_l2;
    
    figure(fig_number)
    plot(freq_ech_X_l2,freq_ech_Y_l2,'d','color',vermillon,'markersize',8,'linewidth',2)    
    nom_figue = [dir_plot 'KIC_' PIC '_echelle_UP.png'];
    set(gcf,'PaperPositionMode','auto')
    eval(['print -dpng -loose ' nom_figue])
   
    %......................................................................
    %                          For final saving
    %......................................................................
    order_n_l0 = round(nu_l0/Dnu - eps_revised);
    order_n_l1 = floor(nu_l1/Dnu - eps_revised - 1/2 + d01_revised); % floor for l=1
    order_n_l2 = round(nu_l2/Dnu - eps_revised - 2/2 + d02);
    
    oscillation.modeID.nu_l0 = nu_l0; oscillation.modeID.n_l0 = order_n_l0; 
    oscillation.modeID.nu_l1 = nu_l1; oscillation.modeID.n_l1 = order_n_l1; 
    oscillation.modeID.nu_l2 = nu_l2; oscillation.modeID.n_l2 = order_n_l2; 
       
    oscillation.M_Msun     = M_Msun;
    oscillation.M_Msun_err = sig_M;
    oscillation.R_Rsun     = R_Rsun;
    oscillation.R_Rsun_err = sig_R;
    
end