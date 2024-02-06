classdef Part < handle 
% Part.m
% Written by Aaron Hess
% University of Oxford
% August 2017
% Data structure for Ensight Gold case file
% Part can only be created from CaseFile.m
% 
% Part class packages data to facililate ensight gold case formated files
% being written
% 
% Use the methods to specify the parts:
%   Geometry
%     - Calculate geometry coordinates for a block of data given a
%       transform (tx). tx is 4x4 affine defining mapping from matlab sub
%       to real space coordinates
%     - connectivity for meshing, ssee connectTypes for allowed tyeps
%     

   properties (Access = private)
      partName;
      
      varList;   % list of variable names
      varHasTime; % boolian list true if var has time steps 
      nTimeStep;
      isGeoBlock;  % is it a structered block (true) or a list of points?
      timeStep;  % contains nTimeStep structure with 
                 % .geometry, 
                 % .connectivity, 
                 % .var is a cell list as long as varList
                 % .szNodes this is 3x1 for isGeoBlock or 1x1 for not isGeoBlock
   end
   
   properties (Constant, Access = public)
      connectTypes = {'point', 'bar2', 'bar3', 'tria3', 'tria6', 'quad4', 'quad8', 'tetra4', 'tetra10', 'pyramid5', 'pyramid13', 'hexa8', 'hexa20', 'penta6', 'penta15'};
   end
   
   methods (Access = {?CaseFile})
       % CaseFile must create a part to ensure all have same n time steps
      function obj = Part(part_name,nstep)
          obj.partName = part_name;
          obj.nTimeStep = nstep;
          obj.timeStep(nstep).geometry = [];
          obj.timeStep(nstep).connectivity = struct();
          obj.timeStep(nstep).var = [];
          obj.timeStep(nstep).szNodes = [];
          obj.isGeoBlock = [];  % not set
          obj.varList = cell(0);
          obj.varHasTime = [];
          
          for iCT = 1:length(obj.connectTypes)              
              obj.timeStep(nstep).connectivity.(obj.connectTypes{iCT}) = [];
          end
      end
   end
   
   methods (Access = public)       
      % Specify coordinates of each node, if block/structered then ndims > 1
      function SetGeometry(obj,tstep,X,Y,Z)
          isBlock = size(X,2)>1;
          if(isempty(obj.isGeoBlock))
             obj.isGeoBlock = isBlock;
          else
              if(isBlock ~= obj.isGeoBlock)
                  if(isBlock)
                      error('SetGeometry Trying to set a structured block geometry when already a non-structered geometry')
                  else
                      error('SetGeometry Trying to set non-structered geometry when already a structered block geometry')
                  end
              end
          end
          if(tstep<1)||(tstep>obj.nTimeStep)
              error('SetGeometry tstep out of range')
          end
          
          if(~isempty(obj.timeStep(tstep).szNodes))
              error('SetGeometry Geometry already set, cannot change')
          end
          
          obj.timeStep(tstep).geometry.X = X;
          obj.timeStep(tstep).geometry.Y = Y;
          obj.timeStep(tstep).geometry.Z = Z;
          
          obj.timeStep(tstep).szNodes = size(X);
      end
      
      % create the coordinates geometry for a unifomrly distributed block
      % according to data size sz
      % tx is the 4x4 affine transform from matlab grid to physical space
      function CreateBlockGeometry(obj,tx,sz)
          if(~isempty(obj.isGeoBlock))
              error('CreateBlockGeometry Geometry allready set')
          end
          [X,Y,Z] = ind2sub(sz,1:prod(sz));  % coordinates for X,Y,Z
          O = ones(size(X));
          vert = tx*[X;Y;Z;O];
          for iT = 1:obj.nTimeStep
              obj.timeStep(iT).geometry.X = vert(1,:);
              obj.timeStep(iT).geometry.Y = vert(2,:);
              obj.timeStep(iT).geometry.Z = vert(3,:);
              obj.timeStep(iT).szNodes = sz;
          end
          obj.isGeoBlock = true;
      end
      
      % does this variable change with time?
      function hasTime = GetVariableHasTime(obj,var_name)
          [~, list_ind] = UpdateList(obj.varList, var_name);
          if(list_ind>length(obj.varList))
              error(['GetVariableHasTime Variable ',var_name,' not found']);
          end
          hasTime = obj.varHasTime(list_ind);
          
      end
      
      % does this variable change with time?
      function hasTime = GetVariableIndHasTime(obj,var_ind)
          hasTime = obj.varHasTime(var_ind);
      end
      
      % all times set the variable (all times must have same number nodes)
      % OR NO Time dimension set for this it must be nNodes (cannot be 3D).
      % var should be either nNodes x nT x 1/2/3 (scalar/2d/vector)
      %                  or nX x nY x nZ x nT x 1/2/3 (scalar/2d/vector)
      function SetVariableTime(obj,var_name, var)
          [out_list, list_ind] = UpdateList(obj.varList, var_name);
          % to test if new variable added, see if list_ind > length(obj.varList)          
          obj.varList = out_list;
          obj.varHasTime(list_ind) = true;
          
          if(ndims(var)>3)
              sz = size(var);
              var = reshape(var,[prod(sz(1:3)),sz(4:end)]);
          else
              sz = size(var);
              var = reshape(var,[prod(sz(1:end)),1]);
          end
              
          sz = size(var);
          nT = obj.nTimeStep;
          if(sz(2)==1)
              obj.varHasTime(list_ind) = false;
              nT = 1;
          elseif(sz(2) ~= obj.nTimeStep)
              error(['Inccorect size of time dimension was ',num2str(sz(2)),' should be ',num2str(obj.nTimeStep)]);
          end
          
          for iT = 1:nT
              if(sz(1) ~= prod(obj.timeStep(iT).szNodes))
                  error(['Inccorect number of nodes ',num2str(sz(1)),' should be ',num2str(prod(obj.timeStep(iT).szNodes)),' for timestep ',num2str(iT)]);
              end
              obj.timeStep(iT).var{list_ind} = squeeze(var(:,iT,:));
          end
      end
      
      % how many variables this part has
      function n = GetTotalVariables(obj)
          n = length(obj.varList);
      end
      
      % test if variable is present
      function isVar = IsVariablePresent(obj,var_name)
          [~, list_ind] = UpdateList(obj.varList, var_name);
          isVar = (list_ind<=length(obj.varList));              
      end
      
      % Get variable by index
      function var = GetVariableByIndex(obj,tstep,ind)
          if(tstep<1)||(tstep>obj.nTimeStep)
              error('GetVariableByIndex Invalid timestep')
          end
          if(ind<1)||(ind>length(obj.varList))
              error('GetVariableByIndex Invalid index')
          end
          var = obj.timeStep(tstep).var{ind};
      end
      
      % extract varaible by name
      function var = GetVariableByName(obj,tstep,var_name)
          [~, list_ind] = UpdateList(obj.varList, var_name);
          if(list_ind>length(obj.varList))
              error(['Variable ', var_name', not found'])
          end
          if(tstep<1)||(tstep>obj.nTimeStep)
              error('Invalid timestep')
          end
          var = obj.timeStep(tstep).var{list_ind};
      end
      
      % Get list of variable names
      function var_names = GetVariableNames(obj)
          var_names = obj.varList;
      end
      
      % get name of variable at index ind
      function var_name = GetVariableName(obj,ind)
          if(ind<1)||(ind>length(obj.varList))
              error('GetVariableName Invalid index')
          end
          var_name = obj.varList{ind}; 
      end
      
      % get the geometry (coordinates) of the part
      function [X,Y,Z] = GetGeometry(obj,tstep)
          if(tstep<1)||(tstep>obj.nTimeStep)
              error('Invalid timestep')
          end
          if(isfield(obj.timeStep(tstep),'geometry')&&isfield(obj.timeStep(tstep).geometry,'X'))
              X = obj.timeStep(tstep).geometry.X;
              Y = obj.timeStep(tstep).geometry.Y;
              Z = obj.timeStep(tstep).geometry.Z;
          else
              X = [];
              Y = [];
              Z = [];
          end
      end
      
      % how many coordinates are there for this time step
      function sz = GetNodeSize(obj,tstep)
          if(tstep<1)||(tstep>obj.nTimeStep)
              error('Invalid timestep')
          end
          sz = obj.timeStep(tstep).szNodes;
      end
      
      % Types of connections can be
      % point bar2 bar3 tria3 tria6 quad4 quad8 tetra4 tetra10 pyramid5 pyramid13
      % hexa8 hexa20 penta6 penta15
      function SetConnectivity(obj,tstep,connect,type)
          if(tstep<1)||(tstep>obj.nTimeStep)
              error('Invalid timestep')
          end
          if(obj.isGeoBlock)
              error('SetConnectivity cannot set connectivty with block geometry')
          end
          % test that index are in range (convert to zero based)
          max_ind = prod(obj.timeStep(tstep).szNodes);
          is_too_big = sum(connect(:)>max_ind)>0;
          is_too_small = sum(connect(:)<1)>0;
          if(is_too_big||is_too_small)
              error('Index out of range');
          end
          % not testing format here ...?
          if(sum(strcmp(obj.connectTypes,type))==0)
              error('SetConnectivity type unknown');
          end
          % todo, check the size is correct...
          obj.timeStep(tstep).connectivity.(type) = connect;
          
      end
      
      % retreive connectivity by type
      function connect = GetConnectivity(obj,tstep,type)
          if(tstep<1)||(tstep>obj.nTimeStep)
              error('Invalid timestep')
          end
          if(isfield(obj.timeStep(tstep).connectivity,type))
              connect = obj.timeStep(tstep).connectivity.(type);
          else
              connect = [];
          end
      end
      
      function isGeo = GetIsStructeredBlock(obj)
          isGeo = obj.isGeoBlock;
      end
   end

end