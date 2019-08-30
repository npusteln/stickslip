%**************************************************************************
% authors: Nelly Pustelnik and Valerie Vidal                              *
% institution: Laboratoire de Physique de l'ENS de Lyon                   *
% date: October 18 2018                                                   *
% License CeCILL-B                                                        *
%**************************************************************************

function [start,stop,taus,taum] = detect_tstartstop(x1,g,vitesse,kstiff,m,fa,th_noise)
th_event = th_noise*10;
Tg=pi*sqrt(m/kstiff);

Fnorm=x1; % normalized force applied on the slider
vp = -(fa*m*g*1e6/kstiff*(diff(Fnorm))-vitesse); % slider velocity
t = linspace(0,length(Fnorm)-1,length(Fnorm))./fa; % time

%disp(['V=',num2str(vitesse),' microns/s']);
%disp(['Vmean=',num2str(mean(vp)),' microns/s']);

% ------------------------------------
% EVENT DETECTION
% ------------------------------------
% detection if slider velocity vp > threshold
ind = find(vp>=th_event);
k=1;  % nb d'evenements de slips
start(k) = ind(1);
stop=[];
for i=2:length(ind)-1;
    if ind(i+1) - ind(i)>1
        stop(k) = ind(i);
        k=k+1;
        start(k) = ind(i+1);
    end
end

% COMPUTATION OF taum=motion time (creep+sliding)
%                taus=stick time
%                mud=dynamic friction coefficient (on slip phases only)
% all points are above the higher threshold = continuous sliding !
if isempty(stop)==1
    taum=Inf; % continuous sliding !
    taus=0;
    mud=mean(Fnorm);
    % indices pour lesquels le patin glisse
    indm=[1:length(vp)];
    % display message
    disp('!!!!!!!!! CONTINUOUS SLIDING !!!!!!!!! (no tstop on figure)')
    % --> go to program end
    
else % go to lower threshold & analyse
    
    % keep the initial finding for checking
    start0=start;
    stop0=stop;
    
    %-----------------------------------------------------------------
    % temps début et fin de slip donné par le seuil bas (noise level)
    %
    % first start & stop (sliding) ----------------------------------
    k=1;
    if start(1)<stop(1) % le signal commence par une phase stick
        while vp(start(k))>th_noise && start(k)>1
            start(k) = start(k)-1;
        end
        while vp(stop(k))>th_noise  && stop(k)<length(vp)
            stop(k) = stop(k)+1;
            if k+1 <= length(start) % il existe un start ensuite (pb fin sinon)
                if stop(k)>=start(k+1) % on fusionne avec le prochain évènement !
                    start(k+1)=[];
                    if k+1 <= length(stop)
                        stop(k)=stop(k+1);
                        stop(k+1)=[];
                    else
                        stop(k)=length(vp); % à supprimer ensuite car pas réel !
                    end
                end
            end
        end
        if stop(k)==length(vp)
            stop(k)=[]; % suppression
        end
    elseif start(1)>stop(1) % le signal commence par une phase de mouvement (creep/slip)
        dum=1;
        while dum==1 % ***
            while vp(stop(k))>th_noise  && stop(k)<length(vp)
                stop(k) = stop(k)+1;
            end
            % find the right next start !
            if stop(1)>=start(1) % on joint l'évènement suivant
                start(1)=[];
                stop(1)=[];
                % et on reboucle sur ***
            elseif start(1)>stop(1) % malgré la correction précédente
                % alors on cherche le vrai start(1)
                while vp(start(1))>th_noise && start(1)>stop(1)
                    start(k) = start(k)-1;
                end
                if start(1)==stop(1) % combine both events ------------
                    % redefine vector "start" by erasing this value
                    start(1)=[];
                    stop(1)=[];
                else % ok exit while loop and go to next event
                    dum=0;
                end
            end
        end
    end
    %
    % get next sliding events ----------------------------------
    while start(k)<start(end) && stop(k)<stop(end)
        k=k+1;
        if start(1)<stop(1) % le signal commence par une phase stick
            while vp(start(k))>th_noise && start(k)>stop(k-1)
                start(k) = start(k)-1;
            end
            while vp(stop(k))>th_noise  && stop(k)<length(vp)
                stop(k) = stop(k)+1;
                if k+1 <= length(start) % il existe un start ensuite (pb fin sinon)
                    if stop(k)>=start(k+1) % on fusionne avec le prochain évènement !
                        start(k+1)=[];
                        if k+1 <= length(stop)
                            stop(k)=stop(k+1);
                            stop(k+1)=[];
                        else
                            stop(k)=length(vp); % à supprimer ensuite car pas réel !
                        end
                    end
                end
            end
            if stop(k)==length(vp)
                stop(k)=[]; % suppression
            end
        elseif start(1)>stop(1) % le signal commence par une phase de mouvement (creep/slip)
            dum=1;
            while dum==1 % ***
                while vp(stop(k))>th_noise  && stop(k)<length(vp)
                    stop(k) = stop(k)+1;
                end
                % find the right next start !
                if stop(k)>=start(k) % on joint l'évènement suivant
                    start(k)=[];
                    stop(k)=[];
                    % et on reboucle sur ***
                elseif start(k)>stop(k) % malgré la correction précédente
                    % alors on cherche le vrai start(k)
                    while vp(start(k))>th_noise && start(k)>stop(k)
                        start(k) = start(k)-1;
                    end
                    if start(k)==stop(k) % combine both events ------------
                        % redefine vector "start" by erasing this value
                        start(k)=[];
                        stop(k)=[];
                    else % ok exit while loop and go to next event
                        dum=0;
                    end
                end
            end
        end
    end
    
    %-----------------------------------------------------------------
    % Compilation des indices "indm" pour lesquels le patin glisse
    % et des temps à l'arrêt taus et en mouvement taum
    %-----------------------------------------------------------------
    continuous_sliding=isempty(stop);
    
    if continuous_sliding ==1 % NO STICK PHASE !!!
        indm=[1:length(vp)];
        taum=Inf; % temps de glissement infini
        taus=0;
        mud=mean(Fnorm);
        % display message
        disp('!!!!!!!!! CONTINUOUS SLIDING !!!!!!!!! (no tstop on figure)')
        
    else  % THERE IS AT LEAST ONE STICK PHASE
        if start(1)<stop(1) % le signal commence par une phase stick
            indm=[start(1):stop(1)]; % 1st sliding
            if start(1)==1 % le sliding avait commencé avant, ne pas le compter !
                taum=[];
                mud=[];
            else
                taum=t(stop(1))-t(start(1));
                mud=mean(Fnorm(start(1):stop(1)));
            end
            taus=[];
            if length(start)==length(stop)
                for i=2:length(start)
                    indm=[indm start(i):stop(i)];
                    taum=[taum t(stop(i))-t(start(i))];
                    mud=[mud mean(Fnorm(start(i):stop(i)))];
                    if i~=length(start) % on ne compte pas le dernier stick [partiel]
                        taus=[taus t(start(i+1))-t(stop(i))];
                    end
                end
            elseif length(start)>length(stop)
                for i=2:length(start)-1
                    indm=[indm start(i):stop(i)];
                    taum=[taum t(stop(i))-t(start(i))];
                    mud=[mud mean(Fnorm(start(i):stop(i)))];
                    taus=[taus t(start(i+1))-t(stop(i))];
                end
                indm=[indm start(end):length(vp)]; % on récupère le bout où ça glisse à la fin pour le %
                % mais non comptabilisé dans taum car partiel
                taus=[taus t(start(end))-t(stop(end))]; % on compte le dernier stick complet
            end
            
        elseif start(1)>stop(1) % le signal commence par une phase de mouvement (creep/slip)
            indm=[1:stop(1)]; % 1st sliding [non comptabilisé dans taum car partiel]
            if length(start)==length(stop)
                for i=2:length(start)
                    indm=[indm start(i-1):stop(i)];
                    if i~=length(start) % on ne compte pas le dernier stick [partiel]
                        taum=[taum t(stop(i+1))-t(start(i))];
                        mud=[mud mean(Fnorm(start(i):stop(i+1)))];
                    end
                    taus=[taus t(start(i))-t(stop(i))];
                end
                indm=[indm start(end):length(vp)]; % on récupère le bout où ça glisse à la fin pour le %
                % mais non comptabilisé dans taum car partiel
            elseif length(start)<length(stop)
                for i=2:length(start)-1
                    indm=[indm start(i-1):stop(i)];
                    taum=[taum t(stop(i+1))-t(start(i))];
                    mud=[mud mean(Fnorm(start(i):stop(i+1)))];
                    taus=[taus t(start(i))-t(stop(i))];
                end
                taus=[taus t(start(end))-t(stop(end))]; % on compte le dernier stick complet
            end
        end
        
    end
    
    % Si taum est vide (des bouts au début ou à la fin du signal)
    if isempty(taum)==1
        taum=Inf; % continuous sliding !
        mud=mean(Fnorm);
    end
    % Si taus est vide (continuous sliding ou bien des bouts au début et à la
    % fin du signal)
    if isempty(taus)==1
        taus=0;
    end
    
end