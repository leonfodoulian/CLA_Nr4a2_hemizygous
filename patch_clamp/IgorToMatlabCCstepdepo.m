   filename=['exp_CCstepDepo_ch10_',num2str(trace),'.ibw'];
   D = IBWread(filename);
   AllTracesCH10(trace,:)=(D.y)';
   ClampedCurrent(trace,1)=((mean(AllTracesCH10(trace,20000:50000)))-(mean(AllTracesCH10(trace,1:9990))))*10^12; % to have it in pA
   filename=['exp_CCstepDepo_ch4_',num2str(trace),'.ibw'];
   D = IBWread(filename);
   FullTraceCH4(trace,:)=(D.y)';
   AllTracesCH4(trace,:)=(D.y(10000:60000))'; %only taking the rec from 0.2 to 1.2 sec 
   [pks,locs,widths,proms]=findpeaks(AllTracesCH4(trace,:),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',200,'MinPeakProminence',0.005,'Annotate','extents');
   nbrAP(trace,:)=[trace,ClampedCurrent(trace,1),length(pks)];
   %locationAP(trace,:)= (locs);
   if ~isempty(widths)  
   ParametersofPeaks(trace,:)=[trace,ClampedCurrent(trace,1),((widths(1)*200)/10000),((mean(widths)*200)/10000),proms(1),mean(proms)]; %widths 1st peak ; mean widths all peaks; amplitude 1st peak; Mean Amp all peaks
   else 
   ParametersofPeaks(trace,:)=[trace,ClampedCurrent(trace,1),0,0,0,0];
   end 
   if ClampedCurrent(trace,1)<-15 && ClampedCurrent(trace,1)>-25
       nm=nm+1;
       HyperpoMinus20pA(nm,:) = FullTraceCH4(trace,:);
       TracesNumberof20pA=trace;
       %plot(HyperpoMinus20pA(nm,:))
   end
   
              if nbrAP(trace,2)<2 
              il=il+1;
              APat0pA(il,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<30 && nbrAP(trace,2)>15
              ij=ij+1;
              APat25pA(ij,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<60 && nbrAP(trace,2)>30
              ik=ik+1;
              APat50pA(ik,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<80 && nbrAP(trace,2)>65
              im=im+1;
              APat75pA(im,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<105 && nbrAP(trace,2)>90
              in=in+1;
              APat100pA(in,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<130 && nbrAP(trace,2)>115
              io=io+1;
              APat125pA(io,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<160 && nbrAP(trace,2)>140
              ip=ip+1;
              APat150pA(ip,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<180 && nbrAP(trace,2)>165
              iq=iq+1;
              APat175pA(iq,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<210 && nbrAP(trace,2)>190
              ir=ir+1;
              APat200pA(ir,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<230 && nbrAP(trace,2)>215
              is=is+1;
              APat225pA(is,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<260 && nbrAP(trace,2)>240
              it=it+1;
              APat250pA(it,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<280 && nbrAP(trace,2)>265
              iu=iu+1;
              APat275pA(iu,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<310 && nbrAP(trace,2)>290
              iv=iv+1;
              APat300pA(iv,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<330 && nbrAP(trace,2)>315
              iw=iw+1;
              APat325pA(iw,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<360 && nbrAP(trace,2)>340
              ix=ix+1;
              APat350pA(ix,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<380 && nbrAP(trace,2)>365
              iy=iy+1;
              APat375pA(iy,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<410 && nbrAP(trace,2)>390
              iz=iz+1;
              APat400pA(iz,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<430 && nbrAP(trace,2)>415
              ia=ia+1;
              APat425pA(ia,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<460 && nbrAP(trace,2)>440
              ib=ib+1;
              APat450pA(ib,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<480 && nbrAP(trace,2)>465
              ic=ic+1;
              APat475pA(ic,:)=nbrAP(trace,3);
              else if nbrAP(trace,2)<510 && nbrAP(trace,2)>490
              id=id+1;
              APat500pA(id,:)=nbrAP(trace,3);
                   
               end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end         
       end
   end   
