% All parameters From seimens protocol (text secion found in dicom header)
% if the field does not exist, then it = 0
%PosSCT:  sSliceArray.asSlice[0].sPosition.dSag
%         sSliceArray.asSlice[0].sPosition.dCor
%         sSliceArray.asSlice[0].sPosition.dTra
%NormSCT: sSliceArray.asSlice[0].sNormal.dSag
%         sSliceArray.asSlice[0].sNormal.dCor
%         sSliceArray.asSlice[0].sNormal.dTra
% This code is taken from the Siemens file fGSLCalcPRS.cpp
function [phase, read, orienation] = fGSLCalcPRS(NormSCT, inPlaneRoation)

%% First Calculate the orientation
  SAGITTAL = 1; CORONAL = 2; TRANSVERSE = 3; % ATH 02/2018 for consistancy
  
  dAbsSagComp     = abs(NormSCT(1));
  dAbsCorComp     = abs(NormSCT(2));
  dAbsTraComp     = abs(NormSCT(3));
  bAlmEqualSagCor = abs(NormSCT(1) - NormSCT(2)) <= 1.e-6;
  bAlmEqualSagTra = abs(NormSCT(1) - NormSCT(3)) <= 1.e-6; %fGSLAlmEqual(dAbsSagComp, dAbsTraComp);
  bAlmEqualCorTra = abs(NormSCT(2) - NormSCT(3)) <= 1.e-6; %fGSLAlmEqual(dAbsCorComp, dAbsTraComp);

  %---------------------------------------------------------------------------
  % Check all values to determine the slice orientation (sag, cor, tra)
  %---------------------------------------------------------------------------
  if ((bAlmEqualSagCor              &&  bAlmEqualSagTra)             ||...
      (bAlmEqualSagCor              &&  (dAbsSagComp < dAbsTraComp)) ||...
      (bAlmEqualSagTra              &&  (dAbsSagComp > dAbsCorComp)) ||...
      (bAlmEqualCorTra              &&  (dAbsCorComp > dAbsSagComp)) ||...
      ((dAbsSagComp > dAbsCorComp)  &&  (dAbsSagComp < dAbsTraComp)) ||...
      ((dAbsSagComp < dAbsCorComp)  &&  (dAbsCorComp < dAbsTraComp)) ||...
      ((dAbsSagComp < dAbsTraComp)  &&  (dAbsTraComp > dAbsCorComp)) ||...
      ((dAbsCorComp < dAbsTraComp)  &&  (dAbsTraComp > dAbsSagComp)))
  
    %-------------------------------------------------------------------------
    % Mainly transverse...
    %-------------------------------------------------------------------------
    
	  orienation = TRANSVERSE;
  
  else
      if ((bAlmEqualSagCor              &&  (dAbsSagComp > dAbsTraComp)) ||...
           (bAlmEqualSagTra              &&  (dAbsSagComp < dAbsCorComp)) ||...
           ((dAbsSagComp < dAbsCorComp)  &&  (dAbsCorComp > dAbsTraComp)) ||...
           ((dAbsSagComp > dAbsTraComp)  &&  (dAbsSagComp < dAbsCorComp)) ||...
           ((dAbsSagComp < dAbsTraComp)  &&  (dAbsTraComp < dAbsCorComp)))
  
           orienation = CORONAL;
      else
          if ((bAlmEqualCorTra              &&  (dAbsCorComp < dAbsSagComp)) ||...
           ((dAbsSagComp > dAbsCorComp)  &&  (dAbsSagComp > dAbsTraComp)) ||...
           ((dAbsCorComp > dAbsTraComp)  &&  (dAbsCorComp < dAbsSagComp)) ||...
           ((dAbsCorComp < dAbsTraComp)  &&  (dAbsTraComp < dAbsSagComp)))
            %-------------------------------------------------------------------------
            % Mainly sagittal...
            %-------------------------------------------------------------------------
            orienation = SAGITTAL;
          else
              %-------------------------------------------------------------------------
              % Invalid slice orientation...
              %-------------------------------------------------------------------------
              error('\n Slice Orientation invalid\n');
          end
      end
  end
  
  switch (orienation)
    case TRANSVERSE
	    phase(1) = 0;
	    phase(2) = NormSCT(3) * sqrt (1. / (NormSCT(2) * NormSCT(2) + NormSCT(3) * NormSCT(3)));  %dGs[2] * sqrt (1. / (dGs[1] * dGs[1] + dGs[2] * dGs[2]));
	    phase(3) = -NormSCT(2) * sqrt (1. / (NormSCT(2)^2 + NormSCT(3)^2));                       %-dGs[1] * sqrt (1. / (dGs[1] * dGs[1] + dGs[2] * dGs[2]));    
 
    case CORONAL
	    phase(1) = NormSCT(2) * sqrt (1. / (NormSCT(1)^2 + NormSCT(2)^2)) ;                       %dGs[1] * sqrt (1. / (dGs[0] * dGs[0] + dGs[1] * dGs[1]));
	    phase(2) = -NormSCT(1) * sqrt (1. / (NormSCT(1)^2 + NormSCT(2)^2));                       %-dGs[0] * sqrt (1. / (dGs[0] * dGs[0] + dGs[1] * dGs[1]));
	    phase(3) = 0.;

      case SAGITTAL
	    phase(1) = -NormSCT(2) * sqrt (1. / (NormSCT(1)^2 + NormSCT(2)^2));                       %-dGs[1] * sqrt (1. / (dGs[0] * dGs[0] + dGs[1] * dGs[1]));
	    phase(2) = NormSCT(1) * sqrt (1. / (NormSCT(1)^2 + NormSCT(2)^2));                        %dGs[0] * sqrt (1. / (dGs[0] * dGs[0] + dGs[1] * dGs[1]));
	    phase(3) = 0.;

  end

  %*------------------------------------------------------------------------*/
  %*  Calculate GR = GS x GP                                                */
  %*------------------------------------------------------------------------*/
  read(1) = NormSCT(2) * phase(3) - NormSCT(3) * phase(2);  % dGs[1] * dGp[2] - dGs[2] * dGp[1];
  read(2) = NormSCT(3) * phase(1) - NormSCT(1) * phase(3);  % dGs[2] * dGp[0] - dGs[0] * dGp[2];
  read(3) = NormSCT(1) * phase(2) - NormSCT(2) * phase(1);  % dGs[0] * dGp[1] - dGs[1] * dGp[0];


  if (inPlaneRoation ~= 0)
      dPhi = inPlaneRoation;
    %*----------------------------------------------------------------------*/
    %* Rotate around the S axis                                             */
    %*----------------------------------------------------------------------*/
	  phase(1) = cos (dPhi) * phase(1) - sin (dPhi) * read(1); %dGp[0] = cos (dPhi) * dGp[0] - sin (dPhi) * dGr[0];
	  phase(2) = cos (dPhi) * phase(2) - sin (dPhi) * read(2); %dGp[1] = cos (dPhi) * dGp[1] - sin (dPhi) * dGr[1];
	  phase(3) = cos (dPhi) * phase(3) - sin (dPhi) * read(3); %dGp[2] = cos (dPhi) * dGp[2] - sin (dPhi) * dGr[2];

    %*----------------------------------------------------------------------*/
    %* Calculate new GR = GS x GP                                           */
    %*----------------------------------------------------------------------*/
	  read(1) = NormSCT(2) * phase(3) - NormSCT(3) * phase(2); % dGr[0] = dGs[1] * dGp[2] - dGs[2] * dGp[1];
	  read(2) = NormSCT(3) * phase(1) - NormSCT(1) * phase(3); % dGr[1] = dGs[2] * dGp[0] - dGs[0] * dGp[2];
	  read(3) = NormSCT(1) * phase(2) - NormSCT(2) * phase(1); %	dGr[2] = dGs[0] * dGp[1] - dGs[1] * dGp[0];
  end

  
  
end