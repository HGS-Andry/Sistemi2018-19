function [mod,fase,omega]=asbode(num,den,W_axes,M_axes,F_axes,B,R,X,Y);
% ASBODE : traccia il diagramma asintotico di Bode 
% ******************************************************************************
%                 ## SYNTAX ##
%
%  ASBODE(NUM,DEN)
%     Traccia il diagramma di Bode di una funzione di trasferimento
%
%                   NUM     b_m*s^m + b_{m-1}*s^{m-1} + ... + b_1*s + b_0
%          W(s) =  ----- = -----------------------------------------------
%                   DEN     a_n*s^n + a_{n-1}*s^{n-1} + ... + a_1*s + a_0
%
%     e scrive sullo schermo i termini che compongono la fattorizzazione di
%     Bode della W(s).
%
%     Il vettore 
%          NUM = [ b_m  b_{m-1}  ...  b_1  b_0 ]  
%     contiene i coefficienti del polinomio al numeratore mentre il vettore
%          DEN = [ a_n  a_{n-1}  ...  a_1  a_0 ]  
%     contiene i coefficienti del polinomio al denominatore della W(s).
%
%     I diagrammi asintotici dei singoli termini sono tracciati in diversi colori,
%     il diagramma asintotico complessivo e' tracciato con una linea nera spessa, 
%     il diagramma reale complessivo e' tracciato con una linea nera tratteggiata. 
%     Le scale e il vettore delle frequenze sono scelti automaticamente.
%
%  ************************************************************************
%  ASBODE(NUM,DEN,[W1,W2],[M1,M2],[F1,F2],B,R,X,Y)
%
%     - [W1,W2] specifica che la scala delle ascisse in entrambi i diagrammi e' 
%        compresa fra 10^W1 e 10^W2 (default [])
%     - [M1,M2] specifica che la scala delle ordinate nel diagramma del modulo e' 
%        compresa fra M1 e M2 (default [])
%     - [F1,F2] specifica che la scala delle ordinate nel diagramma della fase e' 
%       compresa fra F1 e F2 (default [])
%    
%     Il parametro B specifica se calcolare la banda passante
%     - B = 0 : (default) non calcola la banda passante
%     - B = x : (x>0) calcola la banda passante a -x db e indica la
%               corrispondete pulsazione sul diagramma dei moduli
%               NB: la banda si calcola solo per un filtro passa-basso
%                   entro i valori di pulsazione del diagramma
%     
%     Il parametro R specifica se calcolare i punti di risonanza
%     - R = 0 : (default) non calcola modulo, fase e pulsazione alla 
%                risonanza
%     - R = 1 : calcola modulo, fase e pulsazione alla risonanza
%               (punti di massimo locale del diagramma dei moduli
%               entro i valori di pulsazione del diagramma)
%     - R = 2 : calcola modulo, fase e pulsazione alla risonanza
%               indicando tali valori sul diagramma
%     - R = 3 : calcola modulo, fase e pulsazione alla risonanza
%               indicando tali valori solo sul diagramma dei moduli
%     - R = 4 : calcola modulo, fase e pulsazione alla risonanza
%               indicando tali valori solo sul diagramma delle fasi
%
%     Il parametro X specifica su che finestra tracciare i diagrammi
%     - X = 0 : (default) traccia il diagramma dei moduli e delle fasi
%                   nella stessa figura
%     - X = 1 : traccia il diagramma dei moduli e delle fasi in due figure
%    
%     Il parametro Y specifica che tipo di grafico tracciare
%     - Y = 0 : (default) traccia sia i diagrammi asintotici dei singoli
%               termini sia quelli complessivi (reali e asintotici)
%     - Y = 1 : traccia solo i diagrammi reali complessivi
%     - Y = 2 : traccia solo i diagramma asintotici complessivi
%     - Y = 3 : traccia entrambi i diagrammi complessivi
%               - la curva nera tratteggiata indica i diagrammi reali
%               - la curva blu continua indica i diagrammi asintotici
%    
%  ************************************************************************
%  [MOD,FASE,W]=ASBODE(NUM,DEN)
%  [MOD,FASE,W]=ASBODE(NUM,DEN,[W1,W2],[M1,M2],[F1,F2])
%      Non traccia il diagramma ma salva i valori di modulo, fase e frequenza
%      nelle matrici
%      - MOD:  la matrice del moduli in decibel
%              - la prima colonna contiene il diagramma reale complessivo,
%              - la seconda colonna contiene il diagramma asintotico complessivo,
%              - le altre colonne contengono i diagrammi asintotici dei singoli termini
%      - FASE: la matrice delle fasi in gradi
%              - la prima colonna contiene il diagramma reale complessivo,
%              - la seconda colonna contiene il diagramma asintotico complessivo,
%              - le altre colonne contengono i diagrammi asintotici dei singoli termini
%      - W:    il vettore delle frequenze
%
%     Si veda anche BODE.
%
% Autore: Alessandro Giua. Versione 4.1, Maggio 2019.

%%%% settings
RisRealAsLineWidth = 1; %width linea risultante
OthRealAsLineWidth = 0.5; %width altre linee


ni=nargin;
no=nargout;

if (ni<2),
    disp('Errore: la funzione richiede almeno i due argomenti NUM e DEN')
    return
end

if (ni==5), 
   if (size(W_axes) ~= 2),
      disp('Errore: controllare il terzo argomento [W1,W2] (scala delle frequenze)')
      return
   end
   if (size(M_axes) ~= 2),
      disp('Errore: controllare il quarto argomento [M1,M2] (scala dei moduli)')
      return
   end
   if (size(F_axes) ~= 2),
      disp('Errore: controllare il quinto argomento [F1,F2] (scala delle fasi)')
      return
   end
end
if ((ni>=6) & (B<0)),
    disp('Errore: il sesto argomento B puo'' assumere solo valori non negativi')
    return
end
if ((ni>=7) & (R~=0) & (R~=1) & (R~=2) & (R~=3) & (R~=4)),
    disp('Errore: il settimo argomento R puo'' assumere solo valori 0, 1, 2, 3 o 4')
    return
end
if ((ni>=8) & (X~=0) & (X~=1)),
    disp('Errore: l''ottavo argomento X puo'' assumere solo valore 0 e 1')
    return
end
if ((ni>=9) & (Y~=0) & (Y~=1) & (Y~=2) & (Y~=3)),
    disp('Errore: il nono argomento Y puo'' assumere solo valori 0, 1, 2 o 3')
    return
end

if (ni <3),
    W_axes = [];
end
if (ni <4),
    M_axes = [];
end
if (ni <5),
    F_axes = [];
end
if (ni <6),
    B = 0;
end
if (ni <7),
    R = 0;
end
if (ni <8),
    X = 0;
end
if (ni <9),
    Y = 0;
end

% ELIMINO GLI ZERI INIZIALI DEI VETTORI num E den
while ((num(1) == 0) && (max(size(num))>1))
    num = num(2:max(size(num)));
end
while ((den(1) == 0) && (max(size(den))>1))
    den = den(2:max(size(den)));
end



% CALCOLO ZERI E POLI E I COEFFICIENTI NON NULLI DI GRADO PIU' BASSO
z=roots(num);
p=roots(den);
a0=0;
b0=0;
for i=1:length(num),
  if num(i) ~= 0,
      b0=num(i);
  end
end
for i=1:length(den),
  if den(i) ~= 0,
      a0=den(i);
  end
end

% COSTRUISCO UN VETTORE DELLE FREQUENZE CONPOSTO DA 500 PUNTI
% SPAZIATI LOGARITMICAMENTE
if isempty(W_axes),
    [i,j,nonzero]=find([z' p']);
    if isempty(nonzero),
        low=-1;
        up=1;
    else     
    low=floor(log10(0.1*min(abs(nonzero))));
    up=ceil(log10(10*max(abs(nonzero))));
    end
else
    low=W_axes(1);
    up=W_axes(2);
    if ((low ~= round(low)) | (up ~= round(up)) | (low >= up)),
    disp('Errore: per costruire il vettore di frequenze fra 10^W1 e 10^W2')
    disp('        occorre che W1 e W2 siano interi con W1 < W2 (p.e., [-1,2]).')
    return
    end
end
omega=logspace(low,up,500);

% Costruisco una matrice A con tre colonne 
% e in cui ogni riga e' associata ad un fattore
% della forma di Bode
% [K         0      1]   --> caso 1: guadagno
% [nu        0      2]   --> caso 2: nu = eccedenza poli zeri in s=0
% [tau'      0      3]   --> caso 3: zero reale
% [omega_n'  zeta'  4]   --> caso 4: coppia zeri complessi
% [tau       0      5]   --> caso 5: polo reale
% [omega_n   zeta   6]   --> caso 6: coppia poli complessi
% NB: vi e' sempre una e una sola riga con terzo indice 1,
%     vi e' al massimo una riga con terzo indice 2,
%     vi sono zero o piu' righe con terzo indice {3,4,5,6}
m=length(z);
n=length(p);
nu=0;
A=[b0/a0 0 1];
nu_z=length(find(z==0));   % numero di zeri nell'origine 
nu_p=length(find(p==0));   % numero di poli nell'origine 
nu=nu_p-nu_z;              % eccedenza poli-zeri nell'origine
if nu ~= 0,
    A = [A; nu 0 2];
end
for i=1:m,
    if (z(i)~=0 & isreal(z(i))),
        A = [A; -1/z(i) 0 3];
    elseif imag(z(i))>0,
        A = [A; abs(z(i)) -real(z(i))/abs(z(i)) 4];
    end
end
for i=1:n,
    if (p(i)~=0 & isreal(p(i))),
        A = [A; -1/p(i) 0 5];
    elseif imag(p(i))>0,
        A = [A; abs(p(i)) -real(p(i))/abs(p(i)) 6];
    end
end

[ma,na]=size(A);
A=sortrows(A,3);
% SCRIVO SULLO SCHERMO LA FATTORIZZAZIONE DI BODE
for j=1:ma,
   switch A(j,3),
   case 1
       k=A(j,1);
       k_db=round(20*log10(abs(k)));
       if k>0,
        fprintf('\nGuadagno:        K = %7.3f,   K_db = %g db,     phi = 0 deg',k,k_db)
       else
        fprintf('\nGuadagno:        K = %7.3f,   K_db = %g db,     phi = +/- 180 deg',k,k_db)
       end    
   case 2
       nu=A(j,1);
       if nu >0,
           fprintf('\nPoli in origine: nu = %g',nu)
       else
           fprintf('\nZeri in origine: nu'' = %g',-nu)
       end
   case 3
       tau=A(j,1);
       z=-1/tau;
       fprintf('\nZero reale:      z = %7.3f,   tau = %7.3f,   1/|tau| = %7.3f,   phi = da 0 a %+g deg',z,tau,abs(z),sign(tau)*90)
   case 4 
       omega_n=A(j,1);
       zeta=A(j,2);
       re=-zeta*omega_n;
       im=sqrt(omega_n^2-re^2);
       DeltaM = round(20*log10(2*zeta));
       beta=10^(abs(zeta));
       omega_s=omega_n/beta;
       omega_d=omega_n*beta;
       fprintf('\nZeri complessi:  z,z'' = %7.3f +/- j%7.3f,    omega_n = %7.3f,   zeta = %5.2f',re,im,omega_n,zeta)
       fprintf('\n                 beta = %4.1f,                    omega_s = %7.3f,   omega_d = %7.3f,',beta,omega_s,omega_d)
       if abs(zeta)>0,
        fprintf('\n                 phi = da 0 a %+g deg,          Delta M_db = %+g db',sign(zeta)*180,DeltaM)
       else
        fprintf('\n                 phi = da 0 a +/- 180 deg,       Delta M_db = %+g db',DeltaM)
       end    
   case 5
       tau=A(j,1);
       p=-1/tau;
       fprintf('\nPolo reale:      p = %7.3f,   tau = %7.3f,   1/|tau| = %7.3f,   phi = da 0 a %+g deg',p,tau,abs(p),sign(-tau)*90)
   case 6
       omega_n=A(j,1);
       zeta=A(j,2);
       re=-zeta*omega_n;
       im=sqrt(omega_n^2-re^2);
       DeltaM = round(-20*log10(2*zeta));
       beta=10^(abs(zeta));
       omega_s=omega_n/beta;
       omega_d=omega_n*beta;
       fprintf('\nPoli complessi:  p,p'' = %7.3f +/- j%7.3f,    omega_n = %7.3f,   zeta = %5.2f',re,im,omega_n,zeta)
       fprintf('\n                 beta = %4.1f,                    omega_s = %7.3f,   omega_d = %7.3f,',beta,omega_s,omega_d)
       if abs(zeta)>0,
        fprintf('\n                 phi = da 0 a %+g deg,          Delta M_db = %+g db',sign(-zeta)*180,DeltaM)
       else
        fprintf('\n                 phi = da 0 a +/- 180 deg,       Delta M_db = %+g db',DeltaM)
       end    
   end
end
fprintf('\n')

% CALCOLO I DIAGRAMMI ASINTOTICI DEI SINGOLI TERMINI
for i=1:length(omega),
   w=omega(i);
   for j=1:ma,
       switch A(j,3),
       case 1
           k=A(j,1);
           gain(i,j)=20*log10(abs(k));
           phase(i,j)=-180*angle(k)/pi; %costante
       case 2
           nu=A(j,1);
           gain(i,j)=-nu*20*log10(w);
           phase(i,j)=-nu*90;
       case 3
           tau=A(j,1);
           if w < 1/abs(10*tau),
               gain(i,j)=0;
               phase(i,j)=0;
           end
           if ((w >= 1/abs(10*tau)) & (w < 1/abs(tau))),
               gain(i,j)=0;
               phase(i,j)=sign(tau)*45*(log10(w)-log10(1/abs(10*tau)));
           end
           if ((w >= 1/abs(tau)) & (w < 10/abs(tau))),
               gain(i,j)=20*log10(w*abs(tau));
               phase(i,j)=sign(tau)*45*(log10(w)-log10(1/abs(10*tau)));
           end
           if w >= 10/abs(tau),
               gain(i,j)=20*log10(w*abs(tau));
               phase(i,j)=sign(tau)*90;
           end         
       case 4
           omega_n=A(j,1);
           zeta=A(j,2);
           if w < omega_n*10^(-abs(zeta)),
               gain(i,j)=0;
               phase(i,j)=0;
           end
           if ((w > omega_n*10^(-abs(zeta))) & (w < omega_n)),
               gain(i,j)=0;
               phase(i,j)=(90/zeta)*(log10(w/omega_n)+abs(zeta));
           end
           if  w == omega_n,
               gain(i,j)=0;
               phase(i,j)=90*sign(zeta);
           end
           if ((w > omega_n) & (w <= omega_n*10^(abs(zeta)))),
               gain(i,j)=40*log10(w/omega_n);
               phase(i,j)=(90/zeta)*(log10(w/omega_n)+abs(zeta));
           end
           if w > omega_n*10^(abs(zeta)),
               gain(i,j)=40*log10(w/omega_n);
               if zeta ~= 0,
                   phase(i,j)=sign(zeta)*180;
               else
                   phase(i,j)=180;
               end
           end         
       case 5
           tau=A(j,1);
           if w < 1/abs(10*tau),
               gain(i,j)=0;
               phase(i,j)=0;
           end
           if ((w >= 1/abs(10*tau)) & (w < 1/abs(tau))),
               gain(i,j)=0;
               phase(i,j)=-sign(tau)*45*(log10(w)-log10(1/abs(10*tau)));
           end
           if ((w >= 1/abs(tau)) & (w < 10/abs(tau))),
               gain(i,j)=-20*log10(w*abs(tau));
               phase(i,j)=-sign(tau)*45*(log10(w)-log10(1/abs(10*tau)));
           end
           if w >= 10/abs(tau),
               gain(i,j)=-20*log10(w*abs(tau));
               phase(i,j)=-sign(tau)*90;
           end         
       case 6
           omega_n=A(j,1);
           zeta=A(j,2);
           if w < omega_n*10^(-abs(zeta)),
               gain(i,j)=0;
               phase(i,j)=0;
           end
           if ((w > omega_n*10^(-abs(zeta))) & (w < omega_n)),
               gain(i,j)=0;
               phase(i,j)=-(90/zeta)*(log10(w/omega_n)+abs(zeta));
           end
           if  w == omega_n,
               gain(i,j)=0;
               phase(i,j)=-90*sign(zeta);
           end
           if ((w > omega_n) & (w <= omega_n*10^(abs(zeta)))),
               gain(i,j)=-40*log10(w/omega_n);
               phase(i,j)=-(90/zeta)*(log10(w/omega_n)+abs(zeta));
           end
           if w > omega_n*10^(abs(zeta)),
               gain(i,j)=-40*log10(w/omega_n);
               if zeta ~= 0,
                   phase(i,j)=-sign(zeta)*180;
               else
                   phase(i,j)=180;
               end
           end        
       end
   end
end  
% AGGIUNGO ALLE MATRICI GAIN E PHASE:
% - UNA PRIMA COLONNA CHE CONTIENE IL DIAGRAMMA REALE COMPLESSIVO 
% - UNA SECONDA COLONNA CHE CONTIENE IL DIAGRAMMA ASINTOTICO COMPLESSIVO
%AG [G,P]=bode(num,den,omega);
im=sqrt(-1);
for i=1:500,
    G(i)=0;
    P(i)=0;
    w=omega(i);
    for j=1:ma,
        switch A(j,3),
        case 1
            k=A(j,1);
            G(i)=G(i)+20*log10(abs(k));
            P(i)=P(i)+180*angle(k)/pi;
       case 2
           nu=A(j,1);
           G(i)=G(i)-nu*20*log10(w);
           P(i)=P(i)-nu*90;
       case 3
           tau=A(j,1);
           G(i)=G(i)+20*log10(abs(1+im*tau*w));
           P(i)=P(i)+180*angle(1+im*tau*w)/pi;
       case 4
           omega_n=A(j,1);
           zeta=A(j,2);
           G(i)=G(i)+20*log10(abs(1-w^2/omega_n^2+im*2*zeta*w/omega_n));
           P(i)=P(i)+180*angle(1-w^2/omega_n^2+im*2*zeta*w/omega_n)/pi;
       case 5
           tau=A(j,1);
           G(i)=G(i)-20*log10(abs(1+im*tau*w));
           P(i)=P(i)-180*angle(1+im*tau*w)/pi;    
       case 6
           omega_n=A(j,1);
           zeta=A(j,2);
           G(i)=G(i)-20*log10(abs(1-w^2/omega_n^2+im*2*zeta*w/omega_n));
           P(i)=P(i)-180*angle(1-w^2/omega_n^2+im*2*zeta*w/omega_n)/pi;
       end
    end
end
gain=[G' gain*ones(ma,1) gain];
phase=[P' phase*ones(ma,1) phase];
% VALUTO SE OCCORRE SPOSTARE DI +/-360 IL DIAGRAMMA DELLE FASI REALE
% PERCHE' COINCIDA CON QUELLO ASINTOTICO
for i=1:500,
    if phase(i,1) < phase(i,2)-180,
        phase(i,1) = P(i)+360;
    end
    if phase(i,1) > phase(i,2)+180,
        phase(i,1) = P(i)-360;
    end
end
% TRACCIO I DIAGRAMMI
clf(findobj('type','figure','name',1))
clf(findobj('type','figure','name',2))
close(findobj('type','figure','name',1))
close(findobj('type','figure','name',2))
if ((no==0) & (Y==0)),
    %
    % TRACCIO I DIAGRAMMI ASINTOTICI DEI MODULI PER I SINGOLI TERMINI
    if (X==0),
        subplot(2,1,1)
    else
        figure(1)
    end
    % SE K=1 E VI SONO ALTRI FATTORI IL DIAGRAMMA DEI MODULI DEL GUADAGNO 
    % NON VIENE TRACCIATO PERCHE' VALE SEMPRE ZERO E NON CONTRIBUISCE AL 
    % DIAGRAMMA RISULTANTE DEI MODULI
    if ((k==1) & (ma>1)),   
      semilogx(omega,gain(:,4:ma+2),'LineWidth',OthRealAsLineWidth);
    else
      semilogx(omega,gain(:,3:ma+2),'LineWidth',OthRealAsLineWidth);
    end 
    hold on;
    %
    % TRACCIO IL DIAGRAMMA COMPLESSIVO DEI MODULI REALE E ASINTOTICO
    semilogx(omega,gain(:,1),'k--',omega,gain(:,2),'k','LineWidth',RisRealAsLineWidth);
    %
    hold off;
    if (X==0),
        %title('Diagramma di Bode','FontSize',12,'FontName','Times New Roman')        
    else
        title('Diagramma di Bode (modulo)','FontSize',12,'FontName','Times New Roman')
        xlabel('Pulsazione \omega [rad/s]','FontSize',12,'FontName','Times New Roman')
    end
    ylabel('Modulo M [db]','FontSize',12,'FontName','Times New Roman')
    if isempty(M_axes)
        m_max=max(max(gain))+10;
        m_min=min(min(gain))-10;
    else
        m_min=M_axes(1);
        m_max=M_axes(2);
        if m_min >= m_max,
            disp('Errore: per la scala dei moduli [M1,M2] deve valere M1 < M2.')
            return
        end
    end
    axis([min(omega) max(omega) m_min m_max]);
    set(gca,'FontSize',12,'FontName','Times New Roman')
    grid
    %
    % TRACCIO I DIAGRAMMI ASINTOTICI DELLE FASI PER I SINGOLI TERMINI
    if (X==0),
        subplot(2,1,2)
    else
        figure(2)
    end
    % SE K=1 E VI SONO ALTRI FATTORI IL DIAGRAMMA DELLE FASI DEL GUADAGNO 
    % NON VIENE TRACCIATO PERCHE' VALE SEMPRE ZERO E NON CONTRIBUISCE AL 
    % DIAGRAMMA RISULTANTE DELLE FASI
    if ((k==1) & (ma>1)),   
      semilogx(omega,phase(:,4:ma+2),'LineWidth',OthRealAsLineWidth);
    else
      semilogx(omega,phase(:,3:ma+2),'LineWidth',OthRealAsLineWidth);
    end      
    hold on;
    %
    % TRACCIO IL DIAGRAMMA DELLE FASI REALE E ASINTOTICO
    semilogx(omega,phase(:,1),'k--',omega,phase(:,2),'k','LineWidth',RisRealAsLineWidth);
    %
    hold off;
    set(gca,'YTick',-720:45:720,'FontSize',12,'FontName','Times New Roman')
    grid
    if (X==1),
        title('Diagramma di Bode (fase)','FontSize',12,'FontName','Times New Roman')        
    end
    xlabel('Pulsazione \omega [rad/s]','FontSize',12,'FontName','Times New Roman')
    ylabel('Fase \phi [gradi]','FontSize',12,'FontName','Times New Roman')
    if isempty(F_axes)
        f_max=max(max(phase))+20;
        f_min=min(min(phase))-20;
    else
        f_min=F_axes(1);
        f_max=F_axes(2);
        if f_min >= f_max,
            disp('Errore: per la scala delle fasi [F1,F2] deve valere F1 < F2.')
            return
        end
    end
    axis([min(omega) max(omega) f_min f_max]);
end 
if ((no==0) & (Y~=0)), % TRACCIO SOLO I GRAFICI SPECIFICATI DAL PARAMETRO par
    if (X==0),
        subplot(2,1,1)
    else
        figure(1)
    end
    switch Y
          case {1}
            semilogx(omega,gain(:,1),'k','LineWidth',RisRealAsLineWidth);
            if (X==0),
                title('Diagramma di Bode','FontSize',12,'FontName','Times New Roman')
            else
                xlabel('Pulsazione \omega [rad/s]','FontSize',12,'FontName','Times New Roman')
                title('Diagramma di Bode (modulo)','FontSize',12,'FontName','Times New Roman')
            end  
          case {2}
            semilogx(omega,gain(:,2),'b','LineWidth',RisRealAsLineWidth);
            if (X==0),
                title('Diagramma di Bode asintotico','FontSize',12,'FontName','Times New Roman')
            else
                xlabel('Pulsazione \omega [rad/s]','FontSize',12,'FontName','Times New Roman')
                title('Diagramma di Bode asintotico (modulo)','FontSize',12,'FontName','Times New Roman')
            end  
          case {3}
            semilogx(omega,gain(:,1),'k--',omega,gain(:,2),'b','LineWidth',RisRealAsLineWidth);
            if (X==0),
                title('Diagramma di Bode reale e asintotico','FontSize',12,'FontName','Times New Roman')
            else
                xlabel('Pulsazione \omega [rad/s]','FontSize',12,'FontName','Times New Roman')
                title('Diagramma di Bode asintotico (modulo)','FontSize',12,'FontName','Times New Roman')
            end  
    end
    ylabel('Modulo M [db]','FontSize',12,'FontName','Times New Roman')
    if isempty(M_axes),
        m_max=max(max(gain(:,1:2)))+10;
        m_min=min(min(gain(:,1:2)))-10;
    else 
        m_min=M_axes(1);
        m_max=M_axes(2);
        if m_min >= m_max,
            disp('Errore: per la scala dei moduli [M1,M2] deve valere M1 < M2.')
            return
        end
    end
    axis([min(omega) max(omega) m_min m_max]);
    set(gca,'FontSize',12,'FontName','Times New Roman')
    grid
    if (X==0),
        subplot(2,1,2)
    else
        figure(2)
    end
    switch Y
          case {1}
            semilogx(omega,phase(:,1),'k','LineWidth',RisRealAsLineWidth);
            if (X==1),
                title('Diagramma di Bode (fase)','FontSize',12,'FontName','Times New Roman')
            end
          case {2}
            semilogx(omega,phase(:,2),'b','LineWidth',RisRealAsLineWidth);
            if (X==1),
                title('Diagramma di Bode asintotico (fase)','FontSize',12,'FontName','Times New Roman')
            end
          case {3}
            semilogx(omega,phase(:,1),'k--',omega,phase(:,2),'b','LineWidth',RisRealAsLineWidth);
            if (X==1),
                title('Diagramma di Bode reale e asintotico (fase)','FontSize',12,'FontName','Times New Roman')
            end
    end   
    set(gca,'YTick',-720:45:720,'FontSize',12,'FontName','Times New Roman')
    grid
    xlabel('Pulsazione \omega [rad/s]','FontSize',12,'FontName','Times New Roman')
    ylabel('Fase \phi [gradi]','FontSize',12,'FontName','Times New Roman')
    if isempty(F_axes),
        f_max=max(max(phase(:,1:2)))+20;
        f_min=min(min(phase(:,1:2)))-20;
    else
        f_min=F_axes(1);
        f_max=F_axes(2);
        if f_min >= f_max,
            disp('Errore: per la scala delle fasi [F1,F2] deve valere F1 < F2.')
            return
        end
    end
    axis([min(omega) max(omega) f_min f_max]);
end 
if (no>0),
    mod=gain;
    fase=phase;      
end
%
% BANDA PASSANTE
if (B > 0),
    if (nu>0),
        fprintf('\n*** Impossibile calcolare la banda passante: mod(W(0))=Inf\n')
    else
        for i=1:500,
            b(i)=min(0,max((gain(i:500,1)-k_db+B)));
        end
        ibar = find(b,1);
        if isempty(ibar),
            if ((max(size(num)) == max(size(den))) && (num(1)/den(1)<k_db-B)) || (max(size(num)) < max(size(den))),
                fprintf('\n*** Impossibile calcolare la banda passante a -%g db\n    nel range di frequenze selezionato\n',B)
            else
                fprintf('\n*** Impossibile calcolare la banda passante a -%g db\n',B)
            end
        else
            fprintf('\n*** Banda passante a -%g db',B)
            fprintf('\n    Pulsazione: omega_%g = %5.3f rad/s, \t Banda: B_%g = %5.3f Hz\n',B,omega(ibar),B,omega(ibar)/2/pi)
            % Traccio sul diagramma dei moduli a meno che non si tratti del solo diagramma
            % asintotico           
            if ((no==0) & (Y == 2)),
                fprintf('La banda passante non viene indicata sul diagramma asintotico\n')
            end
            if ((no==0) & (Y ~= 2)),
                if (X==0),
                    subplot(2,1,1)
                else
                    figure(1)
                end
                hold on
                semilogx([omega(1),omega(ibar)],[k_db-B,k_db-B],'r:','LineWidth',1.5)
                semilogx([omega(ibar),omega(ibar)],[m_min,gain(ibar,1)],'r:','LineWidth',1.5)
                T1 = strcat('K_{db}-',num2str(B));
                T2 = strcat('\omega_{',num2str(B),'}');
                text(omega(8),k_db-B,T1,'HorizontalAlignment','left','BackgroundColor',[0.9 0.9 0.9])
                text(omega(ibar),m_min+(m_max-m_min)/20,T2,'HorizontalAlignment','center','BackgroundColor',[0.9 0.9 0.9])
                hold off
            end
        end
    end
end
%
% MODULO ALLA RISONANZA
if (R > 0),
    b=zeros(500);
    for i=6:495,
        if (gain(i,1) > max(gain(i-5:i-1,1),gain(i+1:i+5,1))),
            b(i)=1;
        end
    end
    [ibar] = find(b);
    if isempty(ibar),
        fprintf('\n*** Il diagramma non presenta risonanza\n',B)
    else
        h='';
        for i=1:length(ibar),
            w_r = omega(ibar(i));
            M_r = gain(ibar(i),1);
            phi_r = round(phase(ibar(i),1));
            if (length(ibar)==1),
                fprintf('\n*** Il diagramma presenta un punto di risonanza')
                fprintf('\n    Pulsazione: omega_r = %5.3f rad/s, \t Modulo: M_r = %g db,   \t Fase: phi_r = %g deg',w_r,round(M_r),phi_r)
            else
                h= strcat(h,'''');
                if (i==1),
                    fprintf('\n*** Il diagramma presenta %g punti di risonanza',length(ibar))
                end
                fprintf('\n%g - Pulsazione: omega%s_r = %5.3f rad/s, \t Modulo: M%s_r = %g db,   \t Fase: phi%s_r = %g deg',i,h,w_r,h,round(M_r),h,phi_r)
            end
            % Traccio sul diagramma dei moduli a meno che non si tratti del solo diagramma
            % asintotico
           if ((no==0) & (Y == 2) & (R > 1) & (i == length(ibar))),
                fprintf('I parametri alla risonanza non vengono indicati sul diagramma asintotico\n')
            end
            if ((no==0) & (Y ~= 2) & (R > 1)),
                if ((R == 2) | (R == 3)),
                    if (X==0),
                        subplot(2,1,1)
                    else
                        figure(1)
                    end
                    hold on
                    semilogx([omega(1),w_r],[M_r,M_r],'r:','LineWidth',1.5)
                    semilogx([w_r,w_r],[m_min,M_r],'r:','LineWidth',1.5)
                    T1 = strcat('M',h,'_r');
                    T2 = strcat('\omega',h,'_r');
                    text(omega(8),M_r,T1,'HorizontalAlignment','left','BackgroundColor',[0.9 0.9 0.9])
                    text(w_r,m_min+(m_max-m_min)/20,T2,'HorizontalAlignment','center','BackgroundColor',[0.9 0.9 0.9])
                    hold off
                end
                if ((R == 2) | (R == 4)),
                    if (X==0),
                        subplot(2,1,2)
                    else
                        figure(2)
                    end
                    hold on
                    semilogx([omega(1),w_r],[phi_r,phi_r],'r:','LineWidth',1.5)
                    semilogx([w_r,w_r],[f_min,phi_r],'r:','LineWidth',1.5)
                    T1 = strcat('\phi',h,'_r');
                    T2 = strcat('\omega',h,'_r');
                    text(omega(8),phi_r,T1,'HorizontalAlignment','left','BackgroundColor',[0.9 0.9 0.9])
                    text(w_r,f_min+(f_max-f_min)/20,T2,'HorizontalAlignment','center','BackgroundColor',[0.9 0.9 0.9])
                    hold off
                end
            end
        end
        fprintf('\n')
    end
end
fprintf('\n')
return
