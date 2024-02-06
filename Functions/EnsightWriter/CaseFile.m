classdef CaseFile < handle
% CaseFile.m
% Written by Aaron Hess
% University of Oxford
% August 2017
% Encapsulate data to be written into parts
% Usage: c = CaseFile(tSteps) % where tSteps is an array of times for each
%           step
%   Create a new part: p = c.CreatePart('part_name')  % no spaces allowed in name
%   Setup the geometry, connectivity, and variables of the part...
%   c.WriteFile(fname) writes Case files in ASCII format
%

   properties (Access = public)
      partList;  % list of parts, indexed by id
      partNames;
      timeSteps;
      version = '1.0';
   end
   
   methods
      function obj = CaseFile(tSteps)

          if(nargin<1)
              obj.timeSteps = 0;
          else
              a = unique(tSteps);
              if(length(a)<length(tSteps))
                  error('Time Steps must be unique and cannot repeat');
              end
              obj.timeSteps = tSteps;
          end
         
         
         obj.partList=cell(0);
         obj.partNames = cell(0);
      end
      
      % part_name
      function new_part = CreatePart(obj,part_name)
          [namelist, ind] = UpdateList(obj.partNames,part_name);
          
          new_part = Part(part_name,length(obj.timeSteps));
          obj.partList{ind} = new_part;
          obj.partNames = namelist;
          
      end
      
      function part = GetPart(obj,part_name)
          [~, ind] = UpdateList(obj.partNames,part_name);
          if(ind > length(obj.partNames))
              error(['Part ',part_name,' not found']);
          end
          part = obj.partList{ind};
      end
      
      function WriteFile(obj,fname,app)
          [fdir,fbase,~] = fileparts(fname);
          fdir = [fdir,filesep];
          
          if(isempty(obj.partList))
              error('No parts created in case file')
          end
          
          nP = length(obj.partList);
          
          % Write main case file
          fid = fopen(fname,'w');
          if(fid==-1)
              error(['Could not open file ', fname]);
          end
          fprintf(fid,'# Aaron Hess matlab writer version %s\n',obj.version);
          fprintf(fid,'# %s\n',datestr(datetime));
          fprintf(fid,'# Ensight gold case file\n\n');
          
          fprintf(fid,'FORMAT\n');
          fprintf(fid,'type:                      ensight gold\n\n');
          
          fprintf(fid,'GEOMETRY\n');
          fprintf(fid,'model:                     1 %s.geo****\n\n',fbase);
          
          fprintf(fid,'VARIABLE\n');
          var_list = cell(0);
          lastNVar = 0;
          for iP = 1:nP  % for each part print its variables
              for iV = 1:obj.partList{iP}.GetTotalVariables()
                  [var_list,nvar] = UpdateList(var_list,obj.partList{iP}.GetVariableName(iV));
                  if(nvar>lastNVar)  % new variable found, what type is it
                      % record this one exists and what type it is
                      % types can be  scalar, vector, tensor, (not supported complex)
                      % scope can be per case, per node, per element, 
                      sz = size(obj.partList{iP}.GetVariableByIndex(1,iV));
                      switch(sz(2))
                          case 1
                              type = 'scalar';
                          case 3
                              type = 'vector';
                          otherwise
                              type = 'unknown';
                      end
                      if(obj.partList{iP}.GetVariableIndHasTime(iV))
                          fprintf(fid,'%s per node:           1 %s %s.%s****\n',type,var_list{nvar},fbase,var_list{nvar});
                      else
                          fprintf(fid,'%s per node:           1 %s %s.%s0000\n',type,var_list{nvar},fbase,var_list{nvar});
                      end
                          
                  else
                      warning('variables with the same name in different parts not tested');
                  end
                  lastNVar = nvar;
              end
          end
          
          fprintf(fid,'TIME\n');
          fprintf(fid,'time set:               1\n');
          fprintf(fid,'number of steps:        %g\n',length(obj.timeSteps));
          fprintf(fid,'filename start number:  0\n');
          fprintf(fid,'filename increment:     1\n');
          fprintf(fid,'time values:\n');
          for iT = 1:length(obj.timeSteps)
              fprintf(fid,'%e\n',obj.timeSteps(iT));
          end
          
          fid = fclose(fid);
          
          % for each timeStep, write the geometry files and variable files
          %wb = waitbar(0,'Writing Case file');
          progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Please Wait','Message','Writing CASE file'); pause(0.1)

          for iT = 1:length(obj.timeSteps)
              %open geometry file              
              fnameOpen = [fdir,fbase,sprintf('.geo%04g',iT-1)];
              fid = fopen(fnameOpen,'w');
              if(fid == -1)
                  error(['Unable to open file ',fnameOpen]);
              end
              fprintf(fid,'EnSight Model Geometry File\n');
              fprintf(fid,'EnSight 10.2.2\n');
              fprintf(fid,'node id assign\n');
              fprintf(fid,'element id assign\n');
              % calculate extents
              bound = obj.GetGeometryExtents(iT);
              fprintf(fid,'extents\n');
              fprintf(fid,'%12.5e',bound(:,1)); fprintf(fid,'\n');              
              fprintf(fid,'%12.5e',bound(:,2)); fprintf(fid,'\n');              
              fprintf(fid,'%12.5e',bound(:,3)); fprintf(fid,'\n');
              
              for iP = 1:nP
                  fprintf(fid,'part\n');
                  fprintf(fid,'%10d\n',iP);
                  fprintf(fid,'%s\n',obj.partNames{iP});
                  sz=obj.partList{iP}.GetNodeSize(iT);
                  if(obj.partList{iP}.GetIsStructeredBlock())
                      fprintf(fid,'block\n');
                      fprintf(fid,'%10d',sz); fprintf(fid,'\n');
                  else
                      fprintf(fid,'coordinates\n');                      
                      fprintf(fid,'%10d\n',sz(1)); 
                  end
                  [X,Y,Z] = obj.partList{iP}.GetGeometry(iT);
                  fprintf(fid,'%12.5e\n',X(:));
                  fprintf(fid,'%12.5e\n',Y(:));
                  fprintf(fid,'%12.5e\n',Z(:));
                  
                  % Write connectivity
                  cT = obj.partList{iP}.connectTypes;
                  for iCT = 1:length(cT)
                      connect = obj.partList{iP}.GetConnectivity(iT,cT{iCT});
                      if(~isempty(connect))
                          sz = size(connect);
                          fprintf(fid,'%s\n',cT{iCT});
                          fprintf(fid,'%10d\n',sz(1));
                          for iCE = 1:sz(1)
                              fprintf(fid,'%10d',connect(iCE,:)); fprintf(fid,'\n');
                          end
                      end
                  end
              end
              
              fid = fclose(fid);
              
              % for each variable
              for iV = 1:length(var_list)
                  fid = 0;
                  for iP = 1:nP
                      if ((obj.partList{iP}.IsVariablePresent(var_list{iV})) ...
                              && ((iT==1)||obj.partList{iP}.GetVariableHasTime(var_list{iV})))
                          if(fid==0)  %only create file if there is a variable
                              fnameOpen = [fdir,fbase,sprintf('.%s%04g',var_list{iV},iT-1)];
                              fid = fopen(fnameOpen,'w');
                              if(fid == -1)
                                  error(['Unable to open file ',fnameOpen]);
                              end
                              fprintf(fid,'%s\n',var_list{iV});
                          end
                          fprintf(fid,'part\n');
                          fprintf(fid,'%10d\n',iP);
                          if(obj.partList{iP}.GetIsStructeredBlock())
                              fprintf(fid,'block\n');
                          else
                              fprintf(fid,'coordinates\n');
                          end
                          var = obj.partList{iP}.GetVariableByName(iT,var_list{iV});
                          fprintf(fid,'%12.5e\n',var(:));
                      end
                  end
                  if(fid~=0)
                      fid = fclose(fid);
                  end
              end
              % real number are 12.5e format
              % integers are 10d
              %waitbar(iT/length(obj.timeSteps),wb);
              progbar.Value = iT/length(obj.timeSteps);
          end
          %close(wb);
          close(progbar);
      end
      
      function bound = GetGeometryExtents(obj,tstep)
          if(tstep<1)||(tstep>length(obj.timeSteps))
              error('GetGeometryExtents time step out of range');
          end
              
          mins = [1 1 1]*1e10;
          maxs = -mins;
          for iP = 1:length(obj.partList)              
              [X,Y,Z] = obj.partList{iP}.GetGeometry(tstep);
              mi_x = min(X(:));
              ma_x = max(X(:));
              mi_y = min(Y(:));
              ma_y = max(Y(:));
              mi_z = min(Z(:));
              ma_z = max(Z(:));
              if(mins(1)>mi_x)
                  mins(1) = mi_x;
              end
              if(maxs(1)<ma_x)
                  maxs(1) = ma_x;
              end
              if(mins(2)>mi_y)
                  mins(2) = mi_y;
              end
              if(maxs(2)<ma_y)
                  maxs(2) = ma_y;
              end
              if(mins(3)>mi_z)
                  mins(3) = mi_z;
              end
              if(maxs(3)<ma_z)
                  maxs(3) = ma_z;
              end
          end
          bound = [mins;maxs];
      end
   end
%    methods (Access = {?MyOtherClass})
%       function d = myMethod(obj)
%          d = obj.Prop1 + x;
%       end
%    end

end