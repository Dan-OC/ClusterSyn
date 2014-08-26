%%Synthetic signal tests

%another test

%Physical constants:
mi=1.67e-27;
me=9.11e-31;
unit_charge=1.6e-19;
epsilon_0=8.85e-12;
mu_0=4*pi*1e-07;
speedoflight=3e8;

%Plasma parameters
meanB=5e-9;
density=5e5;
omega_ci=(unit_charge*meanB)/mi;
omega_ce=unit_charge*meanB/me;
omega_pi=sqrt((density*(unit_charge)^2)/(mi*epsilon_0));
omega_pe=sqrt((density*(unit_charge)^2)/(me*epsilon_0));
alfven_i=meanB/((mu_0*density*mi)^0.5);

%Times
maxtime=3000;
deltt=0.05;
Ntime=floor(maxtime/deltt);
time=linspace(0,maxtime,Ntime);

%Frequencies:
Nfreq=50;
maxf=1.3*omega_ci/(2*pi);
frequencies=linspace(0.05*omega_ci/(2*pi),maxf,Nfreq);
omega=frequencies*2*pi;

%Probes
chi=1.2e5;
x1=chi;
x2=x1+chi;

%Beall parameters
T=2000;
n_real=10000;
Nk=400;
Nomegalin=400;

%Amplitudes
amp=ones(1,Nfreq); 


%Dispersion relation for left wave:
for i=1:Nfreq
    kioncycl(1,i)=((omega(1,i))/speedoflight)*sqrt(1+(omega_pe^2+omega_pi^2)/( (omega_ce+omega(1,i))*((omega_ci)-((omega(1,i))))));
end

%Signal synthesis
for i=1:Nfreq
signal1(i,1:Ntime)=amp(1,i).*cos(kioncycl(1,i)*x1-omega(1,i).*time);
signal2(i,1:Ntime)=amp(1,i).*cos(kioncycl(1,i)*x2-omega(1,i).*time);
end


%Sum signals
signal1=sum(signal1,1);
signal2=sum(signal2,1);

%Fourier transform:
N=length(signal1);
Y1=fft(signal1)/N;
Y2=fft(signal2)/N;
Fs=1/deltt;
f=Fs/2*linspace(0,1,N/2+1);


%Plot time series:
figure(1);
set(gcf,'Position',[50 50 900 900]);
subplot(3,2,1);
hold off
plot(time, real(signal1),'b',time, real(signal2),'r');
xlabel('Time (s)')
ylabel('Field Amplitude')
xlim([0 T]);
title('Time Series')
legend('C1','C2')
pause(1)
%Plot transform
figure(1)
hold on
subplot(3,2,2);
plot(f*2*pi/omega_ci,2*abs(Y1(1,1:N/2+1)),'b',f*2*pi/omega_ci,2*abs(Y2(1,1:N/2+1)),'r');
xlim([0 1.1*max(omega/omega_ci)]); % 0 1.1*2*max(abs(Y2))])
xlabel('\omega/\omega_{ci}');
ylabel('Amplitude');
title('Periodogram');
legend('C1','C2')
hold off
pause(1)
%Dispersion relations:
figure(1);
subplot(3,2,[3 4]);
hold on
plot(real(kioncycl*alfven_i/omega_ci),omega/omega_ci,'b')
% plot(kioncycr*alfven_i/omega_ci,omega/omega_ci,'g')
% plot(kalf*alfven_i/omega_ci,omega/omega_ci,'m')
xlabel('k V_{A}/\omega_{ci}')
ylabel('\omega/\omega_{ci}')
title('Dispersion Relation')
xlim([0 ceil(max(kioncycl)*alfven_i/omega_ci*1.1)])
ylim([ 0 max(omega)/omega_ci])
hold off
pause(2);



%%Beall analysis

%Split signal into realizations
start=linspace(0,maxtime,n_real);
maxi=find(start<=maxtime-T,1,'last'); %Maximum start point of realization
lengthreal=find(time==start(1,1),1,'first')+ find(time>=T,1,'first') -2;

for i=1:maxi
    
    Tstart=start(i); %Get start point for realization
    Tend=Tstart+T; %get end point
    TstartIdx=find(time>=Tstart,1,'first'); %Get times for these points
    TendIdx=find(time>=Tend,1,'first');
    
    Vf9=signal1(TstartIdx:TendIdx); %Find the signal for these times
    Vf10=signal2(TstartIdx:TendIdx);
       
    if length(Vf9)~=lengthreal %Removes excess if there is any
        Vf9(:,lengthreal+1:lengthreal+(length(Vf9)-lengthreal))=[];
        Vf10(:,lengthreal+1:lengthreal+(length(Vf10)-lengthreal))=[];
    end
    
    NV=length(Vf9); %Calculate how many points, for fft (and defines # of frequencies)
    FsV=1/(time(3)-time(2)); %Nyquist frequency
    fV=FsV/2*linspace(0,1,NV/2+1); %Define the frequencies for realization fft
    omega_sc(i,1:NV/2+1)=2*pi*fV;

    
    fftVf9=fft(Vf9)/NV; %realization fft
    fftVf10=fft(Vf10)/NV;
    
    H(i,:)=conj(fftVf9).*fftVf10;
    S1(i,:)=abs(fftVf9).^2;
    S2(i,:)=abs(fftVf10).^2;
   
    theta=atan2(imag(H(i,:)),real(H(i,:)));
   
    K(i,:)=theta/chi;
   
end


M=size(K,1); %How many realizations there are.
S_hat=zeros(Nomegalin,Nk); %Preallocation
Klin=linspace(0,max(real(kioncycl)),Nk); %Confines K to range 0 to pi/dx CHANGED
omegalin=linspace(0,maxf*2*pi,Nomegalin);
deltaK=Klin(2)-Klin(1); %Difference between K bins
deltaomega=omegalin(2)-omegalin(1);
deltaomegalin=omegalin(2)-omegalin(1);
Nomega=NV/2+1;


% Bin realisations k  (eq 36 from Beall paper)
disp('Binning...')
for OM=1:Nomegalin
 for F=1:Nk
    for omegabin=1:Nomega
       for z=1:M
            if (K(z,omegabin)>=Klin(F) && K(z,omegabin)<Klin(F)+deltaK) && (omega_sc(z,omegabin)>=omegalin(OM) && omega_sc(z,omegabin)<omegalin(OM)+deltaomega)
                %Checks for positive difference of K from mean
                S_hat(OM,F) = S_hat(OM,F) + 0.5*(S1(z,omegabin)+S2(z,omegabin)); %If so, the averages of S1 and S2 are added to S_hat           
            end
       end
    end
 end
end

S_hat=S_hat/M; %Averages out S_hat for number of rows


%Power plot
figure(1);
subplot(3,2,[5 6]);
h2=pcolor(Klin*(alfven_i/omega_ci),omegalin/omega_ci,log10(S_hat)); %Colour plot, dispersion relation with S_hat
ylim([ 0 max(omega)/omega_ci])
xlim([ 0 ceil(max(kioncycl)*alfven_i/omega_ci*1.1)])
hold on
plot(real(kioncycl*alfven_i/omega_ci),omega/omega_ci,'r')
pause(1)
set(h2,'EdgeColor','none');
shading flat;
xlabel('k V_{A}/\omega_{ci}');
ylabel('\omega/\omega_{ci}');
title('Power Distribution, S(K,\omega)');

    
