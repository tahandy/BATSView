classdef BATSView < handle
   
   %-------------------------------------------------------------------
   %                              File struct
   %-------------------------------------------------------------------
   %    fname: Filename 
   %    step:  Step number from CRASH
   %    time:  Simulation time from CRASH
   %    ndim:  Dimensionality of the data
   %    nzone: Number of datapoints in the plotfile
   %    param: Vector storing parameter values
   %    pmap:  Map taking parameter name to its index in 'param
   %    var:   Array storing variable data (size nzone x nvar)
   %    vmap:  Map taking variable name to its index in 'param'
   %    coord: Array of coordinate data (size nzone x ndim)
    
    properties
        % fileObjs: Cell array containing instances of the 'file' struct.
        % nFiles:   Number of loaded files.
        fileObjs;
        nFiles;
    end
        
        
        
    methods
        %===============================================================
        % BATSView:
        %   Construct BATSView object.
        %      
        %   A single BATSView object is allowed to process multiple files
        %   at once. Therefore, filesIn may either be a single filename, or
        %   a cell array of filenames.  
        %
        %   A consequence of loading multiple files is that most member
        %   functions require a fileNum value. For single file instances,
        %   this should be set to 1.
        %
        %   If no input is supplied for variablesIn, all variables are
        %   loaded into the file. As with the filesIn argument, variablesIn
        %   may be either a single variable string or a cell array of
        %   variable strings.
        %
        %   fileTypeIn defaults to the single precision IDL file. The
        %   alternative is 'double', which corresponds to the idl_real8
        %   type in PARAM.in. Ascii has not yet been implemented.
        %
        %===============================================================
        function this = BATSView(filesIn, variablesIn, fileTypeIn)
            
            this.nFiles = 0;
            
            if(nargin>0)
                
                if(iscell(filesIn))
                    files = filesIn;
                elseif(ischar(filesIn))
                    files = {filesIn};
                else
                    error('Unknown case for filesIn\n');
                end
                
                if(~exist('fileTypeIn','var'))
                    fileType = 'single';
                else
                    fileType = fileTypeIn;
                end
                
                if(~exist('variablesIn','var'))
                    variables = 'LOAD_ALL_VAR';
                else
                    if(iscell(variablesIn))
                        variables = variablesIn;
                    elseif(ischar(variablesIn))
                        variables = {variablesIn};
                    else
                        error('[BATSView::constructor] Incapable of handling variables to load\n');
                    end
                end
                
                for i=1:numel(files)
                    this.read_binary(files{i},variables,fileType);
                end
                                
            end
        end
        
        
        %===============================================================
        % get_coords:
        %   Retrieve the coordinate array for a given file number.
        %===============================================================
        function coords = get_coords(this,fileNum)
            if(fileNum > this.nFiles)
                error('[get_coords] File number exceeds number of loaded files\n');
            end
            coords = this.fileObjs{fileNum}.coord;
        end
        
        
        %===============================================================
        % get_var:
        %   Retrieve the variable array for a given file number
        %===============================================================
        function var = get_var(this,fileNum,varname)
            if(fileNum > this.nFiles)
                error('[get_var] File number exceeds number of loaded files\n');
            end
            ind = this.fileObjs{fileNum}.vmap(varname);
            var = this.fileObjs{fileNum}.var(:,ind);
        end
                
        %===============================================================
        % get_uniform:
        %   Provided a subdomain and a minimum spacing size, return a 
        %   uniform node-based representation of the data. Used for
        %   displaying via pcolor.
        %
        %   fileNum: Which loaded file you wish to process.
        %   varname: Name of the variable you wish to process. 
        %   bndbox:  3x2 array identifying the region you wish to map to
        %            the uniform grid. bndbox(:,1) is the lower-left corner
        %            of the box, while bndbox(:,2) is the upper-right
        %            corner.
        %   minsizes: 3x1 array specifying the spacing of the uniform grid
        %            in each direction.
        %===============================================================
        function [X,Y,Z,D] = get_uniform(this,fileNum,varname,bndbox,minsizes)
            if(fileNum > this.nFiles)
                error('[get_var] File number exceeds number of loaded files\n');
            end
            
            ind = this.fileObjs{fileNum}.vmap(varname);
            var = this.fileObjs{fileNum}.var(:,ind);
            
            coord0 = this.get_coords(fileNum);
            ndim   = this.fileObjs{fileNum}.ndim;
            
            nPts = zeros(3,1);
            nPts(1) = ceil((bndbox(1,2)-bndbox(1,1))/minsizes(1));
            if(ndim>1)
                nPts(2) = ceil((bndbox(2,2)-bndbox(2,1))/minsizes(2));
            end
            if(ndim>2)
                nPts(3) = ceil((bndbox(3,2)-bndbox(3,1))/minsizes(3));
            end
            
            switch(ndim)
                case(1)
                    
                    tx = linspace(0,1,nPts(1));
                    xs = bndbox(1,1) + (bndbox(1,2)-bndbox(1,1))*tx;
                    
                    la = true(size(coord0,1),1);
                    la = la & coord0(:,1)>=bndbox(1,1) & coord0(:,1)<=bndbox(1,2);
                     
                    D = interp1(coord0(la,1),var(la),xs,'linear','extrap');
                    X = xs;
                    Y = [];
                    Z = [];
                    
                case(2)
                                        
                    tx = linspace(0,1,nPts(1));
                    ty = linspace(0,1,nPts(2));
                    
                    [TX,TY] = ndgrid(tx,ty);
                    
                    X = bndbox(1,1) + (bndbox(1,2)-bndbox(1,1))*TX;
                    Y = bndbox(2,1) + (bndbox(2,2)-bndbox(2,1))*TY;
                    
                    la = true(size(coord0,1),1);
                    la = la & coord0(:,1)>=bndbox(1,1) & coord0(:,1)<=bndbox(1,2);
                    la = la & coord0(:,2)>=bndbox(2,1) & coord0(:,2)<=bndbox(2,2);
                    
                    SI = scatteredInterpolant(coord0(la,1:2),var(la),'linear','nearest');
                    D = SI(X,Y);
                    Z = [];
                                        
                case(3)
                                        
                    tx = linspace(0,1,nPts(1));
                    ty = linspace(0,1,nPts(2));
                    tz = linspace(0,1,nPts(3));
                    
                    [TX,TY,TZ] = ndgrid(tx,ty,tz);
                    
                    X = bndbox(1,1) + (bndbox(1,2)-bndbox(1,1))*TX;
                    Y = bndbox(2,1) + (bndbox(2,2)-bndbox(2,1))*TY;
                    Z = bndbox(3,1) + (bndbox(3,2)-bndbox(3,1))*TZ;
                    
                    la = true(size(coord0,1),1);
                    la = la & coord0(:,1)>=bndbox(1,1) & coord0(:,1)<=bndbox(1,2);
                    la = la & coord0(:,2)>=bndbox(2,1) & coord0(:,2)<=bndbox(2,2);
                    la = la & coord0(:,3)>=bndbox(3,1) & coord0(:,3)<=bndbox(3,2);
                    
                    SI = scatteredInterpolant(coord0(la,1:3),var(la),'nearest');
                    D = SI(X,Y,Z);
                    
            end
            
        end
        
        %===============================================================
        % get_uniform:
        %   Provided a subdomain and a minimum spacing size, return a 
        %   uniform node-based representation of the data. Used for
        %   displaying via pcolor.
        %
        %   fileNum: Which loaded file you wish to process.
        %   varname: Name of the variable you wish to process. 
        %   bndbox:  3x2 array identifying the region you wish to map to
        %            the uniform grid. bndbox(:,1) is the lower-left corner
        %            of the box, while bndbox(:,2) is the upper-right
        %            corner.
        %   minsizes: 3x1 array specifying the spacing of the uniform grid
        %            in each direction.
        %===============================================================
        function scatint = get_scatteredInterpolant(this,fileNum,varname,bndbox)
            if(fileNum > this.nFiles)
                error('[get_var] File number exceeds number of loaded files\n');
            end
            
            ind = this.fileObjs{fileNum}.vmap(varname);
            var = this.fileObjs{fileNum}.var(:,ind);
            
            coord0 = this.get_coords(fileNum);
            ndim   = this.fileObjs{fileNum}.ndim;
            
%             nPts = zeros(3,1);
%             nPts(1) = ceil((bndbox(1,2)-bndbox(1,1))/minsizes(1));
%             if(ndim>1)
%                 nPts(2) = ceil((bndbox(2,2)-bndbox(2,1))/minsizes(2));
%             end
%             if(ndim>2)
%                 nPts(3) = ceil((bndbox(3,2)-bndbox(3,1))/minsizes(3));
%             end
%             
            switch(ndim)
                case(1)
                    
                    xc      = this.fileObjs{fileNum}.coord(:,1);
                    scatint = griddedInterpolant(xc,var,'pchip');
                    
                case(2)
                                        
%                     tx = linspace(0,1,nPts(1));
%                     ty = linspace(0,1,nPts(2));
%                     
%                     [TX,TY] = ndgrid(tx,ty);
%                     
%                     X = bndbox(1,1) + (bndbox(1,2)-bndbox(1,1))*TX;
%                     Y = bndbox(2,1) + (bndbox(2,2)-bndbox(2,1))*TY;
                    
                    la = true(size(coord0,1),1);
                    la = la & coord0(:,1)>=bndbox(1,1) & coord0(:,1)<=bndbox(1,2);
                    la = la & coord0(:,2)>=bndbox(2,1) & coord0(:,2)<=bndbox(2,2);
                    
                    scatint = scatteredInterpolant(coord0(la,1:2),var(la),'linear','nearest');
                                        
                case(3)
                                        
%                     tx = linspace(0,1,nPts(1));
%                     ty = linspace(0,1,nPts(2));
%                     tz = linspace(0,1,nPts(3));
%                     
%                     [TX,TY,TZ] = ndgrid(tx,ty,tz);
%                     
%                     X = bndbox(1,1) + (bndbox(1,2)-bndbox(1,1))*TX;
%                     Y = bndbox(2,1) + (bndbox(2,2)-bndbox(2,1))*TY;
%                     Z = bndbox(3,1) + (bndbox(3,2)-bndbox(3,1))*TZ;
                    
                    la = true(size(coord0,1),1);
                    la = la & coord0(:,1)>=bndbox(1,1) & coord0(:,1)<=bndbox(1,2);
                    la = la & coord0(:,2)>=bndbox(2,1) & coord0(:,2)<=bndbox(2,2);
                    la = la & coord0(:,3)>=bndbox(3,1) & coord0(:,3)<=bndbox(3,2);
                    
                    scatint = scatteredInterpolant(coord0(la,1:3),var(la),'nearest');                    
            end
            
        end
        
        %===============================================================
        % find_extents:
        %   Given a bounding box, determine the minimum and maximum
        %   (box-like) spatial region, such that value(var)==const. Really
        %   only helpful for integer-like quantities.
        %
        %===============================================================
        function extents = find_extents(this,fileNum,bndbox,varname,varvalue)
            if(fileNum > this.nFiles)
                error('[find_extents] File number exceeds number of loaded files\n');
            end
            
            varind = this.fileObjs{fileNum}.vmap(varname);
            la = abs(this.fileObjs{fileNum}.var(:,varind)-varvalue)<1e-6;
%             la = this.fileObjs{fileNum}.var(:,varind)==0;
            
            
            extents = zeros(3,2);
            
            % x-dim
            labb = (bndbox(1,1)>=this.fileObjs{fileNum}.coord(:,1)) | (this.fileObjs{fileNum}.coord(:,1)>=bndbox(1,2));
            labb = ~labb;
            if(this.fileObjs{fileNum}.ndim>1)
                latmp = (bndbox(2,1)>=this.fileObjs{fileNum}.coord(:,2)) | (this.fileObjs{fileNum}.coord(:,2)>=bndbox(2,2));
                labb = labb & ~latmp;
            end
            if(this.fileObjs{fileNum}.ndim>2)
                latmp = (bndbox(3,1)>=this.fileObjs{fileNum}.coord(:,3)) | (this.fileObjs{fileNum}.coord(:,3)>=bndbox(3,2));
                labb = labb & ~latmp;
            end
                
                
            extents(1,1) = min(this.fileObjs{fileNum}.coord(la & labb,1));
            extents(1,2) = max(this.fileObjs{fileNum}.coord(la & labb,1));
            
            % y-dim
            if(this.fileObjs{fileNum}.ndim>1)
                extents(2,1) = min(this.fileObjs{fileNum}.coord(la & labb,2));
                extents(2,2) = max(this.fileObjs{fileNum}.coord(la & labb,2));
            end
            
            % z-dim
            if(this.fileObjs{fileNum}.ndim>2)
                extents(3,1) = min(this.fileObjs{fileNum}.coord(la & labb,3));
                extents(3,2) = max(this.fileObjs{fileNum}.coord(la & labb,3));
            end
           
        end
            
        
        
        %===============================================================
        % lineout:
        %   Extract a lineout for a given variable. Doesn't do anything
        %   fancy if the line goes outside of the domain, so try to keep it
        %   inside for intelligible (non-extrapolated) results. 
        %
        %   fileNum:  Which loaded file you wish to process.
        %   varname:  Name of the variable you wish to process. 
        %   startIn:  3x1 array specifying the beginning point of the
        %             lineout
        %   finishIn: 3x1 array specifying the ending point of the lineout
        %   nPts:     Number of points to use for the lineout. Uniformly
        %             distributed.
        %
        %===============================================================
        function [coords, values] = lineout(this,fileNum,varname,startIn,finishIn,nPts,bndbox)
            if(fileNum > this.nFiles)
                error('[lineout] File number exceeds number of loaded files\n');
            end
            
            ndim   = this.fileObjs{fileNum}.ndim;
            start  = zeros(3,1);
            finish = zeros(3,1);
            start(1:ndim)  = startIn(1:ndim);
            finish(1:ndim) = finishIn(1:ndim);
            
            ind = this.fileObjs{fileNum}.vmap(varname);
            var = this.fileObjs{fileNum}.var(:,ind);
            
%             bndbox = zeros(3,2);
%             bndbox(1,:) = [min(start(1),finish(1)), max(start(1),finish(1))];
%             bndbox(2,:) = [min(start(2),finish(2)), max(start(2),finish(2))];
%             bndbox(3,:) = [min(start(3),finish(3)), max(start(3),finish(3))];
            
            t = linspace(0,1,nPts)';
            
            switch(ndim)
                case(1)
                    coords = start(1) + (finish(1)-start(1))*t;
                    xc     = this.fileObjs{fileNum}.coord(:,1);
                    values = interp1(xc,var,coords,'nearest');                    
                    
                case(2)
                    coords = zeros(nPts,2);
                    coords(:,1) = start(1) + (finish(1)-start(1))*t;
                    coords(:,2) = start(2) + (finish(2)-start(2))*t;
                    
                    coord0 = this.fileObjs{fileNum}.coord(:,1:ndim);
                    la = true(size(coord0,1),1);
                    if(exist('bndbox','var'))
                        la = la & coord0(:,1)>=bndbox(1,1) & coord0(:,1)<=bndbox(1,2);
                        la = la & coord0(:,2)>=bndbox(2,1) & coord0(:,2)<=bndbox(2,2);
                    end
                                        
                    SI = scatteredInterpolant(this.fileObjs{fileNum}.coord(la,1:2),var(la),'nearest');
                    values = SI(coords);
                    
                case(3)
                    coords = zeros(nPts,3);
                    coords(:,1) = start(1) + (finish(1)-start(1))*t;
                    coords(:,2) = start(2) + (finish(2)-start(2))*t;
                    coords(:,3) = start(3) + (finish(3)-start(3))*t;
                    SI = scatteredInterpolant(this.fileObjs{fileNum}.coord(:,1:3),var,'nearest');
                    values = SI(coords);
                    
            end
            
        end
        
        
        
        %===============================================================
        % lineout_exact:
        %   Extract a lineout for a given variable. Doesn't do anything
        %   fancy if the line goes outside of the domain, so try to keep it
        %   inside for intelligible (non-extrapolated) results. 
        %
        %   fileNum:  Which loaded file you wish to process.
        %   varname:  Name of the variable you wish to process. 
        %   startIn:  3x1 array specifying the beginning point of the
        %             lineout
        %   finishIn: 3x1 array specifying the ending point of the lineout
        %   nPts:     Number of points to use for the lineout. Uniformly
        %             distributed.
        %
        %===============================================================
        function [XY, data] = lineout_exact(this,fileNum,pt,dir,limits,vars)
            
            if(fileNum > this.nFiles)
                error('[lineout] File number exceeds number of loaded files\n');
            end
            
            xmin = limits(1);
            xmax = limits(2);
            y0 = pt(2);

            
            XYZ   = this.get_coords(1);
            x     = XYZ(:,dir);
            dx    = this.get_var(1,'dx');
            
            dxmin = min(dx(:));
            dxmax = max(dx(:));
            
            xCentsNew = [];
            switch(this.fileObjs{1}.ndim)
                case(1)
                    xCentsNew = x;
                case(2)
            
                    bb      = [[xmin,xmax];[y0, y0+dxmax];[-1e99,1e99]];
                    scintU  = this.get_scatteredInterpolant(1,'dx',bb);

                    blksz = 8;
                    x0 = xmin;

                    while(true)

                        lvl_test = log2(dxmax/dxmin);
                        foundBlk = false;
                        for L = lvl_test:-1:0
                            dxloc = dxmin*2^L;
                            x1 = x0 + blksz*dxloc;

                            xf = linspace(x0,x1,blksz+1);
                            xc = 0.5*(xf(2:end)+xf(1:end-1));
                            dxSamp = scintU(xc,y0*ones(size(xc)));

                            if(all(abs(dxSamp-dxSamp(1))<1.0*dxmin) && abs(dxSamp(1)-dxloc)<1.0e0*dxmin)
                                xCentsNew = [xCentsNew, xc];
                                x0 = x1;
                                foundBlk = true;
                                break;
                            end
                        end

                        if(~foundBlk)
                            L
                            [x0, x1]
                            [dxmin*2^L,dxSamp]
                            error('Unable to find block!');
                        end

                        if(abs(x0-xmax)<dxmin)
                            break;
                        end

                    end
                otherwise
                    error('3D not supported yet');
            end

            XY = zeros(numel(xCentsNew),2);
            XY(:,1) = xCentsNew;
            XY(:,2) = pt(2)*ones(size(xCentsNew));
            data = zeros(numel(xCentsNew),numel(vars));
            for i=1:numel(vars)
                switch(this.fileObjs{1}.ndim)
                    case(1)
                        data(:,i) = this.get_var(1,vars{i});
                    case(2)    
                        scint = this.get_scatteredInterpolant(1,vars{i},bb);
                        data(:,i) = scint(XY(:,1),XY(:,2));                
                    otherwise
                        error('3D not supported yet');
                end
            end
            
        end
        
        
        %===============================================================
        % read_binary:
        %   Read a single/double precision plot file produced by BATSRUS.
        %   This may be called repeatedly to load multiple outputs into a
        %   single BATSView object.        
        %===============================================================
        function read_binary(this,filename,variables,prec)
            
            % Set the precision for real-type values
            floatFmt   = 'float32';
            floatBytes = 4;
            if(strcmpi(prec,'double'))
                floatFmt   = 'float64';
                floatBytes = 8;
            end
            
            if(ischar(variables) && strcmpi(variables,'LOAD_ALL_VAR'))
                loadVarSubset = false;
            elseif(iscell(variables))
                loadVarSubset = true;
            else
                error('[read_binary] What did you send as variables?\n');
            end
                               
                
            %-----------------------------------------------------------
            % Begin reading the Fortran binary file. Note that if binary
            % files are written in sequential mode, Fortran adds an
            % additional header and footer to each record. The header is
            % 32-bit integer corresponding to the number of 4-byte words
            % contained in the record. This value is repeated at the end of
            % the record.
            %-----------------------------------------------------------
            
            fprintf('[read_binary] Loading file %s\n',filename);
            
            % Open the provided file
            fid = fopen(filename,'r');
            
            fread(fid,1,'ubit32');
            hdr = fread(fid,500,'ubit8=>char');
            fread(fid,1,'ubit32');

            fread(fid,1,'ubit32');
            step = fread(fid,1,'bit32=>int');
            time = fread(fid,1,floatFmt);
            ndim = fread(fid,1,'bit32=>int');
            npar = fread(fid,1,'bit32=>int');
            nvar = fread(fid,1,'bit32=>int');
            fread(fid,1,'ubit32');

            fread(fid,1,'ubit32');
            n_D = fread(fid,abs(ndim),'bit32=>int');
            fread(fid,1,'ubit32');

            if(npar>0)
                fread(fid,1,'ubit32');
                parv = fread(fid,npar,'float32');
                fread(fid,1,'ubit32');
            end

            fread(fid,1,'ubit32');
            namevar = fread(fid,500,'ubit8=>char')';
            fread(fid,1,'ubit32');

            fread(fid,1,'ubit32');
            coord = fread(fid,[prod(n_D),abs(ndim)],floatFmt);
            fread(fid,1,'ubit32');

            
            
            % Split out the variable names into their own strings.
            % Note that the first 1:nDim variable names are x,[y,[z]],
            % and the last end-npar+1:end names are the parameter names.
            % These are removed from the list of variable names.
            varNames = regexp(strtrim(namevar),'\s+','split');
            varNames(1:abs(ndim)) = [];
            parNames = varNames(end-npar+1:end);
            varNames(end-npar+1:end) = [];
            
            
            if(loadVarSubset)
                vmap = containers.Map();
            else
                vmap = containers.Map(varNames,1:nvar);
            end
            
            nvarout = nvar;
            if(loadVarSubset)
                nvarout = numel(variables);
            end
            
            cnt = 0;
            var = zeros(prod(n_D),nvarout);
            for i=1:nvar                
                fread(fid,1,'ubit32');
                
                if(loadVarSubset)
                    if(any(strcmp(varNames{i},variables)))
                        cnt = cnt + 1;
                        var(:,cnt) = fread(fid,prod(n_D),floatFmt);
                        vmap(varNames{i}) = cnt;
                    else
                        fseek(fid,prod(n_D)*floatBytes,'cof');
                    end
                else
                    cnt = cnt + 1;
                    var(:,cnt) = fread(fid,prod(n_D),floatFmt);
                end
                    
                fread(fid,1,'ubit32');
            end

            fclose(fid);
            
            
                
            fileObj = struct('fname', filename, ...
                             'step',  step, ...
                             'time',  time, ...
                             'ndim',  abs(ndim), ...
                             'nzone', size(var,1), ...
                             'param', parv, ...
                             'pmap',  containers.Map(parNames,1:npar), ...
                             'nvar',  nvarout, ...
                             'var',   var, ...
                             'vmap',  vmap, ...
                             'coord', coord);
                         
            fileObj.x = coord(:,1);
            if(abs(ndim)>1)
                fileObj.y = coord(:,2);
                if(abs(ndim)==3)
                    fileObj.z = coord(:,3);
                end
            end
                
            % Store file object
            this.nFiles = this.nFiles + 1;
            this.fileObjs{this.nFiles} = fileObj;
            
        end
        
        
%         %===============================================================
%         % read_binary:
%         %   Read a single/double precision plot file produced by BATSRUS.
%         %   This may be called repeatedly to load multiple outputs into a
%         %   single BATSView object.        
%         %===============================================================
%         function read_ascii(this,filename,variables,prec)
%             
%             % Set the precision for real-type values
%             floatFmt   = 'float32';
%             floatBytes = 4;
%             if(strcmpi(prec,'double'))
%                 floatFmt   = 'float64';
%                 floatBytes = 8;
%             end
%             
%             if(ischar(variables) && strcmpi(variables,'LOAD_ALL_VAR'))
%                 loadVarSubset = false;
%             elseif(iscell(variables))
%                 loadVarSubset = true;
%             else
%                 error('[read_binary] What did you send as variables?\n');
%             end
%                                
%                 
%             %-----------------------------------------------------------
%             % Begin reading the Fortran binary file. Note that if binary
%             % files are written in sequential mode, Fortran adds an
%             % additional header and footer to each record. The header is
%             % 32-bit integer corresponding to the number of 4-byte words
%             % contained in the record. This value is repeated at the end of
%             % the record.
%             %-----------------------------------------------------------
%             
%             fprintf('[read_ascii] Loading file %s\n',filename);
%             
%             [hdr,data] = mhdrload(filename);
%             var    = containers.Map(regexp(strtrim(hdr(2,:)),'\s+','split'),1:size(data,2));
%             
%             
% TIME   = var('t');
% CHIMAX = var('chimax');
% FOAMIN = var('foammin');
% SHCKMIN = var('shockmin');
% SHCKMAX = var('shockmax');
% DENSL = var('rholeft');
% DENSR = var('rhoright');
% 
% ACCPE = var('accelpex');
% ACCPI = var('accelpix');
%            
        
    end
    
end