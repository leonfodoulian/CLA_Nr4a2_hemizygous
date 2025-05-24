%% Add the .ibw extension to the file missing it

TraceList = dir('exp_*');
for FieldToCheck=1:length(TraceList)
    if endsWith(TraceList(FieldToCheck).name,'ibw')
    else 
    [~, f,ext] = fileparts(TraceList(FieldToCheck).name);
    rename = strcat(f,'.ibw',ext) ; 
    movefile(TraceList(FieldToCheck).name, rename); 
    end 
end 

%% Find the trace number Depo

%clear;clc;

TraceNumberDepo=[];
TraceListCCstepDepo = dir('exp_CCstepDepo_ch4*');
for fieldsDepo=1:length(TraceListCCstepDepo)
    TNDepo=TraceListCCstepDepo(fieldsDepo).name;
    TNstringDepo=convertCharsToStrings(TNDepo(19:24)); %Only works for trace up to #999
    TNCorrDepo=extractBetween(TNstringDepo,"_",".");
    TraceNumberDepo(1,fieldsDepo)=TNCorrDepo;
    TraceNumberDepo=sort(TraceNumberDepo);
end 

TraceNumberHyperpo=[];
TraceListCCstepHyperpo = dir('exp_CCstepHyperpo_ch4*');
for fieldsHyperpo=1:length(TraceListCCstepHyperpo)
    TNHyperpo=TraceListCCstepHyperpo(fieldsHyperpo).name;
    TNstringHyperpo=convertCharsToStrings(TNHyperpo(22:26)); %Only works for trace up to #999
    TNCorrHyperpo=extractBetween(TNstringHyperpo,"_",".");
    TraceNumberHyperpo(1,fieldsHyperpo)=TNCorrHyperpo;
    TraceNumberHyperpo=sort(TraceNumberHyperpo);
end 


TraceNumberAPThreshold=[];
TraceListAPThreshold = dir('exp_APThreshold_ch4*');
for fieldsAPThreshold=1:length(TraceListAPThreshold)
    TNAPThreshold=TraceListAPThreshold(fieldsAPThreshold).name;
    TNstringAPThreshold=convertCharsToStrings(TNAPThreshold(20:24)); %Only works for trace up to #999
    TNCorrAPThreshold=extractBetween(TNstringAPThreshold,"_",".");
    TraceNumberAPThreshold(1,fieldsAPThreshold)=TNCorrAPThreshold;
    TraceNumberAPThreshold=sort(TraceNumberAPThreshold);
end  


%% CCstep depo
%clear;clc;

if ~isempty(TraceListCCstepDepo)

il=0;ij=0;ij=0;ik=0;im=0;in=0;io=0;ip=0;iq=0;ir=0;is=0;it=0;iu=0;iv=0;iw=0;ix=0;
iy=0;iz=0;ia=0;ib=0;ic=0;id=0;ie=0;ig=0; nm=0;


%Enter trace number of depo
warning ('off','all');


prompt="Do you want to Lower peak detection threshold Y/N [Y]:";
txt = input(prompt,"s");
if txt=='Y' | txt=='y'
    prompt="Which threshold would you like to use ?";
    MinPeakHeight=input(prompt);
else 
    MinPeakHeight=0;
end 

prompt="Do you want to use specific traces for CC step DEPO Y/N [Y]:";
txt = input(prompt,"s");
if txt=='Y' | txt=='y'
    prompt="which traces do you want to start with ?";
    tracestart=input(prompt);
    prompt="which traces do you want to end with ?";
    tracestop=input(prompt);


    for trace=tracestart:tracestop
         
    end 
  
else 
   
    MinPeakHeight=0;
    for  T=1:length(TraceNumberDepo)
      trace=[];
      trace=TraceNumberDepo(1,T); 
      IgorToMatlabCCstepdepo
    end
end 

currents = 0:25:500;
MeanFreq_ACSF = zeros(length(currents), 1);

for i = 1:length(currents)
    current = currents(i);
    data = eval(sprintf('APat%dpA', current));
    MeanFreq_ACSF(i) = mean(data);
end

datefilename=pwd;
date= datefilename(end-11:end);
freq='_ACSF_Freq';
filename=strcat(date,freq);
 
save(filename,'MeanFreq_ACSF')
%save('workspaceACSF')


% Find Reobase
kj=1;
if sum(nbrAP(:,3))>0
  while (nbrAP(kj,3))==0
    kj=kj+1;
  end 
   RheoBase=nbrAP(kj,2); 
else 
  RheoBase=NaN; 
end
end  
%% CC step hyperpo 

if ~isempty(TraceListCCstepDepo)
if isempty(TraceListCCstepHyperpo) 
    %if exist('HyperpoMinus20pA') 
 NS= size(HyperpoMinus20pA);
 NbrofSteps=NS(1);
   if NbrofSteps>1 
   MeanHyperpoMinus20pA=mean(HyperpoMinus20pA(:,:));
   else 
   MeanHyperpoMinus20pA=(HyperpoMinus20pA(:,:));
   end 
   figure(1)
   plot(MeanHyperpoMinus20pA) 
else 
   
  prompt="Do you want to use specific traces for CC step HYPERPO Y/N [Y]:";
  txt = input(prompt,"s");
if txt=='Y' | txt=='y'
    prompt="which traces do you want to start with ?";
    tracestart=input(prompt);
    prompt="which traces do you want to end with ?";
    tracestop=input(prompt);
for trace= tracestart:tracestop
    
    filename=['exp_CCstepHyperpo_ch10_',num2str(trace),'.ibw'];
   D = IBWread(filename);
   AllTracesHyperpoCH10(trace,:)=(D.y)';
   ClampedCurrentHyperpo(trace,1)=((mean( AllTracesHyperpoCH10(trace,20000:50000)))-(mean(AllTracesHyperpoCH10(trace,1:9990))))*10^12; % to have it in pA
   filename=['exp_CCstepHyperpo_ch4_',num2str(trace),'.ibw'];
   D = IBWread(filename);
   AllTracesHyperpoCH4(trace,:)=(D.y)'; %only taking the rec from 0.2 to 1.2 sec 
   
   
   if ClampedCurrentHyperpo(trace,1)<-15 && ClampedCurrentHyperpo(trace,1)>-25
       nm=nm+1;
       HyperpoMinus20pA(nm,:) = AllTracesHyperpoCH4(trace,:);
       TracesNumberof20pA=trace;
       %plot(HyperpoMinus20pA(nm,:))
   end
    
end 
else
    
  
   for T=1:length(TraceNumberHyperpo)
   trace=[];
   trace=TraceNumberHyperpo(1,T);
   filename=['exp_CCstepHyperpo_ch10_',num2str(trace),'.ibw'];
   D = IBWread(filename);
   AllTracesHyperpoCH10(trace,:)=(D.y)';
   ClampedCurrentHyperpo(trace,1)=((mean( AllTracesHyperpoCH10(trace,20000:50000)))-(mean(AllTracesHyperpoCH10(trace,1:9990))))*10^12; % to have it in pA
   filename=['exp_CCstepHyperpo_ch4_',num2str(trace),'.ibw'];
   D = IBWread(filename);
   AllTracesHyperpoCH4(trace,:)=(D.y)'; %only taking the rec from 0.2 to 1.2 sec 
   
   
   if ClampedCurrentHyperpo(trace,1)<-15 && ClampedCurrentHyperpo(trace,1)>-25
       nm=nm+1;
       HyperpoMinus20pA(nm,:) = AllTracesHyperpoCH4(trace,:);
       TracesNumberof20pA=trace;
       %plot(HyperpoMinus20pA(nm,:))
   end
end 
end 
end 
end
%end 

 IgortoMatlabCCstepHyperpo
 if exist('HyperpoMinus20pA')
 save('CellIntrinsicProperties','RestingMembranePot', 'DeltaMP','InputResistance','Tauinms','Capacitance','RheoBase','Width1stPeak','MeanWidthAllPeaks','Amp1stPeak','MeanAmpAllPeaks')
 %save('workspace')
 end 
 %savefig(figure(1),'Fig step-20pA With 1 Tau')
 savefig(figure(2),'Fig exponential fitting')
 

%% Timing Prop
 
if ~exist('TraceNumberDepo') 
    TraceNumberDepo=TraceNumber; 
end

if exist('tracestart')
    clear TraceNumberDepo;
    bh=[];
    bh= tracestart:tracestop;
    for count=1:length(bh)
    TraceNumberDepo(1,count)= bh(1,count);
   end 
end
  
TraceNumberDepo=sort(TraceNumberDepo);   %This paragraph doesn't exist in other scripts, let's check if it's really necessary

if TraceNumberDepo(1,1) == (TraceNumberDepo(1,2) - 1)
    disp('CCstep depo 1 does not exist we are good to continue');
else 
    TraceNumberDepo(1, 1:end-1) = TraceNumberDepo(1, 2:end);
    TraceNumberDepo(:, end) = [];
end 


for  T=1:length(TraceNumberDepo)
            trace=[];
      trace=TraceNumberDepo(1,T); 
      %filename=['exp_CCstepDepo_ch4_',num2str(trace),'.ibw'];
       TracePoints(T,:)=FullTraceCH4(trace,[10100:59900]); %to avoid peaks before or after the CC step 
       Derivative(T,:)= diff(TracePoints(T,:));
       y=TracePoints(T,:);
     % Settings
      lag = 30;
      threshold = 4;
      influence = 0.5;
      x1=[1:length(TracePoints)];
 
      [signals,avg,stdFilter] = ThresholdingAlgo(y,lag,threshold,influence,T);

 signals=signals';
 OnOffSignals(T,:)=signals(1,:);  %Save signals result in a matrix
 stdFilter=stdFilter';
 StdDev(T,:)=stdFilter(1,:);
 avg=avg';
 AvgDev(T,:)=avg(1,:);
 
 [pks,locs,widths,proms]=findpeaks(StdDev(T,:),'MinPeakHeight',0,'MinPeakDistance',250,'MinPeakProminence',0.0009,'Annotate','extents');
      maxpeaks{1,T}=locs;
      maxvalue{1,T}=pks;
      amplitude{1,T}=proms;
      
 %[pks,locs,widths,proms]=findpeaks(TracePoints(T,:),'MinPeakHeight',0,'MinPeakDistance',200,'MinPeakProminence',0.005,'Annotate','extents');
 %   MaxPeaksFoundOnRealTrace{1,T}=locs;
 
 
    if isempty(maxpeaks{1,T})
         onsets{1,T}=[];
         offsets{1,T}=[];
    end
   
     
      for k=1:length(maxpeaks{1,T})  
      IndexOnsets= maxpeaks{1,T}(k);
     while StdDev(T,IndexOnsets) > mean(StdDev(T,IndexOnsets-(1:10)))
         IndexOnsets=IndexOnsets-1;
     end 
     while TracePoints(T,IndexOnsets)< mean(TracePoints(T,IndexOnsets+(1:5))) && IndexOnsets < 10
           IndexOnsets=IndexOnsets+1;
      end 
     onsets{1,T}(k) = IndexOnsets;
     end 
      
     for j=1:length(maxpeaks{1,T})
         if maxpeaks{1,T}(j) < 49500
     IndexOffsets= maxpeaks{1,T}(j)+150; %to avoid the double peak 
      while StdDev(T,IndexOffsets) > mean(StdDev(T,IndexOffsets+(1:10))) && IndexOffsets < 49790
     IndexOffsets=IndexOffsets+1;
      end 
      while TracePoints(T,IndexOffsets)> mean(TracePoints(T,IndexOffsets+(1:10))) && IndexOffsets < 49790
           IndexOffsets=IndexOffsets+1;
      end 
       offsets{1,T}(j) = IndexOffsets;
         else
             offsets{1,T}(j)=maxpeaks{1,T}(j);
         end 
     end 
   
     %Let's correct maxpeak who is a bit to early compared to the real peak
   
     for ds=1:length(maxpeaks{1,T})
        windowtoFindMP=[];
        if maxpeaks{1,T}(ds) < 49500
        windowtoFindMP=[maxpeaks{1,T}(ds):maxpeaks{1,T}(ds)+100];
        [CMV,CMI]=max(TracePoints(T,windowtoFindMP)); %Max value and its indice
        CorrectedMaxValue{1,T}(ds)=CMV;
        CorrectedMaxpeakIndice{1,T}(ds)=maxpeaks{1,T}(ds)+(CMI-1); % changed with the CMI-1 use to be just CMI
        else 
         CorrectedMaxValue{1,T}(ds)=TracePoints(T,maxpeaks{1,T}(ds));
         CorrectedMaxpeakIndice{1,T}(ds)=maxpeaks{1,T}(ds);
        end 
    end 
 

if ~isempty(onsets{1,T}); 
    
    for l=1:length(CorrectedMaxpeakIndice{1,T})
    RiseTime{1,T}(l) = ((CorrectedMaxpeakIndice{1,T}(l)) - (onsets{1,T}(l)))/50; %in ms
    if RiseTime{1,T}(l)<0 | RiseTime{1,T}(l)>5 % to remove aberant RT hence all the rest variables
       RiseTimes{1,T}(l)=NaN;
       APTotalTime{1,T}(l)=NaN; 
       AmpRT{1,T}(l)=NaN;
       AmpDT{1,T}(l)=NaN;
       SpeedRT{1,T}(l)=NaN;
       SpeedDT{1,T}(l)=NaN;
    else   
    DecayTime{1,T}(l) = (offsets{1,T}(l)- CorrectedMaxpeakIndice{1,T}(l))/50;
    APTotalTime{1,T}(l)= RiseTime{1,T}(l)+DecayTime{1,T}(l);
    AmpRT{1,T}(l) = (abs(TracePoints(T,onsets{1,T}(l))-TracePoints(T,CorrectedMaxpeakIndice{1,T}(l))))*1000;
    AmpDT{1,T}(l) = (abs(TracePoints(T,offsets{1,T}(l))-TracePoints(T,CorrectedMaxpeakIndice{1,T}(l))))*1000;
    SpeedRT{1,T}(l) = AmpRT{1,T}(l)/RiseTime{1,T}(l);
    SpeedDT{1,T}(l) = AmpDT{1,T}(l)/DecayTime{1,T}(l);
    end
    
   
    %Find the halfwidth time
    
if ~isnan(AmpRT{1,T}(l))
        HalfAmpValue{1,T}(l)= TracePoints(T,onsets{1,T}(l)) + ((AmpRT{1,T}(l)))/2000;
        %HalfAmpValueDrugs{1,T}(l)=((AmpRTDrugs{1,T}(l)))/2000;
        windowtofindminonset =[];
        windowtofindminonset=[onsets{1,T}(l):CorrectedMaxpeakIndice{1,T}(l)];
        windowsizeonset=size(windowtofindminonset);
        Diffonset=[];
        for counter=1:windowsizeonset(2)
          %Diffonset(counter)=abs(abs(TracePoints(T,windowtofindminonset(counter)))-((HalfAmpValue{1,T}(l))));
          Diffonset(counter)= abs(HalfAmpValue{1,T}(l)-(TracePoints(T,windowtofindminonset(counter))));
       
        end 
          [MinvalueOnset,IndexHalfAmponset]=min(Diffonset);
          IndexHalfAmpOnset=(IndexHalfAmponset-1)+onsets{1,T}(l);
          IndexHalfAmpOnsetReal{1,T}(l)=IndexHalfAmpOnset;
          
        windowtofindminoffset =[];
        windowtofindminoffset=[CorrectedMaxpeakIndice{1,T}(l):offsets{1,T}(l)];
        windowsizeoffset=size(windowtofindminoffset);
        for counterbis=1:windowsizeoffset(2)
          %Diffoffset(counter)=(TracePoints(T,windowtofindminoffset(counter)))-(HalfAmpValue{1,T}(l));
          Diffoffset(counterbis)= abs(HalfAmpValue{1,T}(l)-(TracePoints(T,windowtofindminoffset(counterbis))));
        end 
          [MinvalueOffset,IndexHalfAmpoffset]=min(Diffoffset);
          IndexHalfAmpOffset=(IndexHalfAmpoffset-1)+CorrectedMaxpeakIndice{1,T}(l);
          IndexHalfAmpOffsetReal{1,T}(l)=IndexHalfAmpOffset;
          
          HalfWidthTime{1,T}(l)=((IndexHalfAmpOffsetReal{1,T}(l))-(IndexHalfAmpOnsetReal{1,T}(l)))/50;
    else 
         IndexHalfAmpOnsetReal{1,T}(l)=NaN;
         IndexHalfAmpOffsetReal{1,T}(l)=NaN;
         HalfWidthTime{1,T}(l)=NaN;
     
    
    
    % Find offset as y=onsets
    
    
     %TracePoints(T,onsets{1,T}(l))
     windowtofindshortoffset=[];
     windowtofindshortoffset=[CorrectedMaxpeakIndice{1,T}(l):offsets{1,T}(l)];
      wfso=size(windowtofindshortoffset);
        for counterbisbis=1:wfso(2)
          %Diffoffset(counter)=(TracePoints(T,windowtofindminoffset(counter)))-(HalfAmpValue{1,T}(l));
          DiffShortOffset(counterbisbis)= abs(TracePoints(T,onsets{1,T}(l))-(TracePoints(T,windowtofindshortoffset(counterbisbis))));
        end 
        [MinvalueShortOffset,IndexShortoffset]=min(DiffShortOffset);
          IndexShortOffset=(IndexShortoffset-1)+CorrectedMaxpeakIndice{1,T}(l);
          if IndexShortOffset>49801
              IndexShortOffset=49801;
          end 
          
          IndexShortOffsetReal{1,T}(l)=IndexShortOffset;
          MinvalueShortOffsetReal{1,T}(l)=TracePoints(T,IndexShortOffset);
          MinvalueOffsetReal{1,T}(l) =TracePoints(T,offsets{1,T}(l));
          % Amp AHP
          AmpAHP{1,T}(l)=MinvalueShortOffsetReal{1,T}(l)-MinvalueOffsetReal{1,T}(l); %in mV 
          TimeAHP{1,T}(l)=(offsets{1,T}(l)-IndexShortOffsetReal{1,T}(l))/50; %in ms
          SpeedAHP{1,T}(l)=AmpAHP{1,T}(l)/TimeAHP{1,T}(l); % in mV/ms
          
    end
    
    end 
        
end

   if isempty(maxpeaks{1,T}) 
       RiseTime{1,T} = NaN;
       DecayTime{1,T}=NaN;
       APTotalTime{1,T}=NaN; 
       AmpRT{1,T}=NaN;
       AmpDT{1,T}=NaN;
       SpeedRT{1,T}=NaN;
       SpeedDT{1,T}=NaN;
       HalfWidthTime{1,T}=NaN;
       AmpAHP{1,T}=NaN;
       TimeAHP{1,T}=NaN;
       SpeedAHP{1,T}=NaN;
          
   end 
end 
 
for  T=1:length(CorrectedMaxpeakIndice)

for l=1:length(CorrectedMaxpeakIndice{1,T})


    windowtofindshortoffset=[];
     windowtofindshortoffset=[CorrectedMaxpeakIndice{1,T}(l):offsets{1,T}(l)];
      wfso=size(windowtofindshortoffset);
      DiffShortOffset=[];
        for counterbisbis=1:wfso(2)
          %Diffoffset(counter)=(TracePoints(T,windowtofindminoffset(counter)))-(HalfAmpValue{1,T}(l));
          DiffShortOffset(counterbisbis)= abs(TracePoints(T,onsets{1,T}(l))-(TracePoints(T,windowtofindshortoffset(counterbisbis))));
        end 
        [MinvalueShortOffset,IndexShortoffset]=min(DiffShortOffset);
          IndexShortOffset=(IndexShortoffset-1)+CorrectedMaxpeakIndice{1,T}(l);
          IndexShortOffsetReal{1,T}(l)=IndexShortOffset;
          MinvalueShortOffsetReal{1,T}(l)=TracePoints(T,IndexShortOffset);
          MinvalueOffsetReal{1,T}(l) =TracePoints(T,offsets{1,T}(l));
          % Amp AHP
          AmpAHP{1,T}(l)=(MinvalueShortOffsetReal{1,T}(l)-MinvalueOffsetReal{1,T}(l))*1000; %in mV 
          TimeAHP{1,T}(l)=(offsets{1,T}(l)-IndexShortOffsetReal{1,T}(l))/50; %in ms
          SpeedAHP{1,T}(l)=AmpAHP{1,T}(l)/TimeAHP{1,T}(l); % in mV/ms
          AmpDecayShortOffset{1,T}(l) =  (abs(MinvalueShortOffsetReal{1,T}(l))+TracePoints(T,CorrectedMaxpeakIndice{1,T}(l)))*1000;
          DecayTimeShortOffset{1,T}(l) = (IndexShortOffsetReal{1,T}(l) - CorrectedMaxpeakIndice{1,T}(l))/50;
          SpeedDecayTimeShortOffset{1,T}(l) = AmpDecayShortOffset{1,T}(l)/DecayTimeShortOffset{1,T}(l); % in mV/ms
        if AmpAHP{1,T}(l)<0 || AmpAHP{1,T}(l)> 15 || DecayTimeShortOffset{1,T}(l)==0
            AmpAHP{1,T}(l)=[];
            TimeAHP{1,T}(l)=[];
            SpeedAHP{1,T}(l)=[];
            AmpDecayShortOffset{1,T}(l)=[];
            DecayTimeShortOffset{1,T}(l)=[];
            SpeedDecayTimeShortOffset{1,T}(l)=[];
        end 
end 

end 


MeanAmpAHP=nanmean(cell2mat(AmpAHP));
MeanTimeAHP=nanmean(cell2mat(TimeAHP));
MeanSpeedAHP=nanmean(cell2mat(SpeedAHP));
MeanAmpDecaySO=nanmean(cell2mat(AmpDecayShortOffset));
MeanDecayTimeSO=nanmean(cell2mat(DecayTimeShortOffset));
MeanSpeedDecayTimeSO=nanmean(cell2mat(SpeedDecayTimeShortOffset));
MeanRiseTime=nanmean(cell2mat(RiseTime));
MeanDecayTime=nanmean(cell2mat(DecayTime));
MeanAPTotalTime= nanmean(cell2mat(APTotalTime));
MeanAmpRT= nanmean(cell2mat(AmpRT));
MeanAmpDT= nanmean(cell2mat(AmpDT));
MeanSpeedRT=nanmean(cell2mat(SpeedRT));
MeanSpeedDT=nanmean(cell2mat(SpeedDT));
MeanHalfWidthTime=nanmean(cell2mat(HalfWidthTime));

%% Rheobase 

   rd=1;
   SizeParameters=size(ParametersofPeaks);
   IndexColumn=3;
     idxnon0 = find(ParametersofPeaks ~= 0, 1); 
 ParametersofPeaksDepo=ParametersofPeaks(idxnon0:end,:);
 PPD=ParametersofPeaksDepo(:,3)
 % Find indices where the value changes from 0 to non-zero
indicesRehobase = find(PPD(1:end-1) == 0 & PPD(2:end) ~= 0) + 1;
for kp=1:length(indicesRehobase)
Rheobaseformean(kp)=ParametersofPeaksDepo(indicesRehobase(kp),2)
end
 Rheobase=mean(nonzeros(Rheobaseformean))
 
if exist('Rheobaseformean')
Rheobase=mean(nonzeros(Rheobaseformean))
  save('Rheobase','Rheobase')
end
   
%%

  np=0;
for  T=1:length(ClampedCurrent)
      trace=[];
      trace=T; 
      
  if ClampedCurrent(trace,1)<40 && ClampedCurrent(trace,1)>10
       np=np+1;
       Depo25pA(np,:) = FullTraceCH4(trace,:);
       TracesNumberof25pA=trace;
       %plot(HyperpoMinus20pA(nm,:))
  end
end 


if exist('Depo25pA')
   NS= size(Depo25pA);
   NbrofSteps=NS(1);

   if NbrofSteps>1 
   MeanDepo25pA=mean(Depo25pA(:,:));
   else 
   MeanDepo25pA=(Depo25pA(:,:));
   end 
   %plot(MeanHyperpoMinus20pADrugs)
   
   RestingMembranePot25pA=mean(MeanDepo25pA(1:9000))*1000;
   SteadyStateMembranePot25pA=mean(MeanDepo25pA(25000:55000))*1000;
   DeltaMP25pA= abs(SteadyStateMembranePot25pA-RestingMembranePot25pA); 
   InputResistance25pA=(DeltaMP25pA/20)*1000; %in Megaohm
 PercentofMax25pA=(1-exp(-1))*DeltaMP25pA ;
  PotentialOf1Tau25pA=(RestingMembranePot25pA-PercentofMax25pA)/1000;
  
  MinMatrix25pA(1,:)=abs(MeanDepo25pA(10000:15000)-PotentialOf1Tau25pA);
  [M,I] = min(MinMatrix25pA); 
  IndexOfTau25pA=10000+I;
  
  %plot(IndexOfTauDrugs,PotentialOf1TauDrugs,'*r')
  
  Tauinms25pA=(I*200)/10000;  %in ms (10000 point are 200ms)
  Capacitance25pA=Tauinms25pA/InputResistance25pA;   %in nF
        
  Curvetofit25pA=MeanDepo25pA(10000:20000)';
  Curvetofit25pA=-Curvetofit25pA; 
  Curvetofitat025pA=Curvetofit25pA-min(Curvetofit25pA);
  
  TimeMatrix25pA(1,1)=0.02;
  for k=2:length(Curvetofitat025pA)
      TimeMatrix25pA(1,k)=TimeMatrix25pA(1,(k-1))+0.02;
  end
   
  TimeMatrix25pA=TimeMatrix25pA';
  
  [ExpFit,gof3]=fit(TimeMatrix25pA,Curvetofitat025pA,'exp1');
  if gof3.rsquare<0.8
      disp('CAREFULL rsquare is below 80%')
      carefulsentence=('CAREFULL rsquare is below 80%');
      save('BADEXPFITat25pA','carefulsentence')
  end 
  
  Tauinms25pA=-1/(ExpFit.b);
  Capacitance25pA=Tauinms25pA/InputResistance25pA; %in nF
  
  figure(5)
   plot(ExpFit,TimeMatrix25pA,Curvetofitat025pA)
   RS=num2str(gof3.rsquare);
   legend('rsquare=',RS)
   title('Exp fit ACSF 25pA protocol')
end
  
   if exist('Depo25pA')
 save('CellIntrinsicProperties25pA','RestingMembranePot25pA', 'DeltaMP25pA','InputResistance25pA','Tauinms25pA','Capacitance25pA')
 savefig(figure(5),'Fig Exp Fit at 25pA')
   end 

    %% ISI
 
     % Define the range and step size
startValue = 0;
endValue = 500;
stepSize = 25;

% Create the column vector
 MeanFreq_ACSF(:,2) = linspace(startValue, endValue, (endValue - startValue) / stepSize + 1)';
 MinPeakHeight=0; 
  % count=1;
   % while MeanFreq_ACSF(count,1)<7 && count<21
  %     count=count+1;
  %   end
     
    % if count<21
         
     %FirstCCsteptolookforISI_ACSF=MeanFreq_ACSF(count,2)
     FirstCCsteptolookforISI_ACSF=Rheobase*0.7+Rheobase
     if FirstCCsteptolookforISI_ACSF>490
         FirstCCsteptolookforISI_ACSF=500;
     end 
     
     ClampedCurrentdepo = ClampedCurrent(idxnon0:end); 
     tracenumbertofindISI_ACSF= find(ClampedCurrentdepo>FirstCCsteptolookforISI_ACSF-15 & ClampedCurrentdepo<FirstCCsteptolookforISI_ACSF+15);
     
     for hp=1:length(tracenumbertofindISI_ACSF)
         traceISI=tracenumbertofindISI_ACSF(hp,1);
         [pks,locs,widths,proms]=findpeaks(AllTracesCH4(traceISI,:),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',200,'MinPeakProminence',0.005,'Annotate','extents');
         if length(locs)>6
         ISI1_ACSF(hp)=(locs(1,2)-locs(1,1))*200/10000; %in ms
         ISI2_ACSF(hp)=(locs(1,3)-locs(1,2))*200/10000; %in ms
         ISI3_ACSF(hp)=(locs(1,end-1)-locs(1,end-2))*200/10000;
         ISI4_ACSF(hp)=(locs(1,end)-locs(1,end-1))*200/10000; %in ms
         AmpPeak1_ACSF=AmpRT{1,traceISI}(1);
         AmpPeak2_ACSF=AmpRT{1,traceISI}(2);
         AmpPeak3_ACSF=AmpRT{1,traceISI}(3);
         MeanAmpFirst3Peaks=(AmpPeak1_ACSF+AmpPeak2_ACSF+AmpPeak3_ACSF)/3
         AmpPeak4_ACSF=AmpRT{1,traceISI}(end-2);
         AmpPeak5_ACSF=AmpRT{1,traceISI}(end-1);
         AmpPeak6_ACSF=AmpRT{1,traceISI}(end);
         MeanAmpLast3Peaks=(AmpPeak4_ACSF+AmpPeak5_ACSF+AmpPeak6_ACSF)/3
         end 
     end 
     
     if exist ('ISI1_ACSF')
     MeanISI1_ACSF=mean(ISI1_ACSF);
     MeanISI2_ACSF=mean(ISI2_ACSF);
     MeanISI3_ACSF=mean(ISI3_ACSF);
     MeanISI4_ACSF=mean(ISI4_ACSF);
     
     ISIAverage_ACSF=(MeanISI1_ACSF+ MeanISI2_ACSF+MeanISI3_ACSF+MeanISI4_ACSF)/4;
     ISI1RatioISIAverage_ACSF=MeanISI1_ACSF/ISIAverage_ACSF;
     AdaptationIndex_ACSF=((MeanISI1_ACSF+MeanISI2_ACSF)/2)/((MeanISI3_ACSF+MeanISI4_ACSF)/2);
     
     AccomodationIndex_ACSF=MeanAmpFirst3Peaks/MeanAmpLast3Peaks;
     
     save('ISI_ACSF_70%Rheobase','AdaptationIndex_ACSF','AccomodationIndex_ACSF','MeanISI1_ACSF','MeanISI2_ACSF','MeanISI3_ACSF','MeanISI4_ACSF','ISIAverage_ACSF','ISI1RatioISIAverage_ACSF')
     %save workspaceISI_ACSF         
     else 
         EmptyVariables=[0,0,0,0,0,0,0];
         save('ISI_ACSF_70%Rheobase','EmptyVariables')
     end 
     
%% DAP & Delay 1st AP

% Find the rheobase trace 
% Make sure that if the cell is bursting at rheobase we cannot try to find
% a DAP as we won't find one. 

for jh=1:length(indicesRehobase)
    
    offsets{1,indicesRehobase(jh)}
    onsetDAP(jh)=offsets{1,indicesRehobase(jh)}(1);  
    [minDAP(jh),IndminDAP(jh)]=min(TracePoints(indicesRehobase(jh), onsetDAP(jh)-50:onsetDAP(jh)+50));
    IndminDAP(jh)=onsetDAP(jh)-(51-IndminDAP(jh));
    [maxDAP(jh),IndmaxDAP(jh)]= max(TracePoints(indicesRehobase(jh), onsetDAP(jh):onsetDAP(jh)+300));
    IndmaxDAP(jh)=onsetDAP(jh)+IndmaxDAP(jh)-1;
    AmpDAP(jh)=(maxDAP(jh)-minDAP(jh))*1000; %in mV
    onset1stAP(jh)= onsets{1,indicesRehobase(jh)}(1);
    Delay1stAP(jh)= (onset1stAP(jh)/50); %in ms
    
    if length(maxpeaks{1,indicesRehobase(jh)})>1
        if maxpeaks{1,indicesRehobase(jh)}(2)-maxpeaks{1,indicesRehobase(jh)}(1)<5500
            AmpDAP(jh)=0;
        end 
    end 
    
end 
  MeanAmpDAP=mean(nonzeros(AmpDAP));
  MeanDelay1stAP=mean(Delay1stAP);
  
  save('DAP&Delay1stAP','MeanAmpDAP','MeanDelay1stAP')
  
%% Sigmoid Fit
     
FreqCell=MeanFreq_ACSF(:,1);
InjectedCurrent=[0:25:500]';

Maxfreq=max(FreqCell);

% Define the sigmoid model
sigmoidModel = @(params, x) params(1) ./ (1 + exp(-params(2) * (x - params(3))));

% Provide multiple sets of initial guesses
initialGuesses = [
    [1, 1, 1];
    [0.5, 0.5, 0.5];
    [Maxfreq, 0, 200];
    [Maxfreq, 0, 250];
    % Add more sets of initial guesses as needed
];

% Fit the data with multiple initial guesses
bestParams = [];
bestResiduals = Inf;

for i = 1:size(initialGuesses, 1)
    [params, residuals] = lsqcurvefit(sigmoidModel, initialGuesses(i, :), InjectedCurrent, FreqCell); % Use y_double here

    % Check if this fit has smaller residuals
    if norm(residuals) < norm(bestResiduals)
        bestParams = params;
        bestResiduals = residuals;
    end
end


% Define the sigmoid model
sigmoidModel = @(A, k, x0, x) A ./ (1 + exp(-k * (x - x0)));


% Initial guesses
initialGuess = bestParams;

fitresult = fit(InjectedCurrent, FreqCell, sigmoidModel, 'StartPoint', initialGuess);

% Get the residuals
residuals = FreqCell - feval(fitresult, InjectedCurrent);

% Calculate SSR and SST
SSR = sum(residuals.^2);
mean_y = mean(FreqCell);
SST = sum((FreqCell - mean_y).^2);

% Calculate R-squared
R2 = 1 - SSR / SST;


% Visualize the best fit
figure(2)
plot(InjectedCurrent, FreqCell, 'b.');
hold on;
plot(fitresult, 'r-');
legend('Data', 'Fitted Sigmoid Curve');
xlabel('Injected Current (pA)');
ylabel('Frequency (Hz)');
annotation('textbox', [0.15, 0.7, 0.2, 0.2], 'String', ['R-squared: ', num2str(R2)], 'EdgeColor', 'none', 'FontSize', 12, 'Color', 'k');

% Display the best fit parameters
disp(['Best Fit Parameters: ', num2str(bestParams)]);

if R2<0.95
    BADSIGMOIDFIT_ACSF='BADSIGMOIDFIT_ACSF'
    save(BADSIGMOIDFIT_ACSF,'BADSIGMOIDFIT_ACSF')
end 

%save('ParamsSigmoidDrugs.mat','bestParamsDrugs')
%savefig(figure(2),'FigSigmoidDrugs')

save('R2SigmoidACSF.mat','R2')

%% SAG and rebound

ClampedCurrentHyperpo(ClampedCurrentHyperpo == 0) = NaN;
ClampedCurrentHyperpo_Rounded=round(ClampedCurrentHyperpo/10)*10;
ClampedCurrentHyperpo_Rounded_NonNaN = ClampedCurrentHyperpo_Rounded(~isnan(ClampedCurrentHyperpo_Rounded));
uniqueGroupsHyperpo = unique(ClampedCurrentHyperpo_Rounded_NonNaN); 

indicesHyperpoSteps = cell(length(uniqueGroupsHyperpo), 1);

for i = 1:length(uniqueGroupsHyperpo)
    indicesHyperpoSteps{i} = find(ClampedCurrentHyperpo_Rounded == uniqueGroupsHyperpo(i));
end

MeanPerStepHyperpo = cell(length(uniqueGroupsHyperpo), 1);

for i = 1:length(uniqueGroupsHyperpo)
    if length(indicesHyperpoSteps{i,1})>1
    MeanPerStepHyperpo{i,1} = mean(AllTracesHyperpoCH4(indicesHyperpoSteps{i,1},:));
    else 
    MeanPerStepHyperpo{i,1} = AllTracesHyperpoCH4(indicesHyperpoSteps{i,1},:);
    end 
end

MinSAG=min(MeanPerStepHyperpo{1,1}(10000:15000)); %find exact min just after the hyperpo step of -120pA
MaxSAG=mean(MeanPerStepHyperpo{1,1}(45000:59000)); %the max of the sag is the mean of the 15000
SAG=(MaxSAG-MinSAG)*1000; %to have it in mV

MaxRebound=max(MeanPerStepHyperpo{1,1}(62000:75000)); 
MinRebound=min(MeanPerStepHyperpo{1,1}(65000:75000)); 
Rebound=(MaxRebound-MinRebound)*1000; 

save('SAG&Rebound','SAG','Rebound')
