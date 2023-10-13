classdef GPT_FacetInjector_NMM < handle & physConsts & batchLSF
  % GPT_FACETINJECTOR
  
  properties

    %NMM: Uninteresting setup stuff that cannot be deleted
    tmpDirectory='/scratch';
    jobName='F2Injector';
    dataDirectory=pwd;%'/home/glen/OneDrive/Work/FACET2/data/F2Injector';
    runDirectory=pwd;%'/home/glen/OneDrive/Work/FACET2/data/F2Injector';% CE change 6/2023
    batchQueueName='short'; % nmacro<=1e5 (short) 1e5<nmacro<1e7 (medium)
    mcrRoot='/nfs/slac/g/nlcrd/u01/whitegr/mcr/v96';
    srcdir=pwd%'/home/glen/OneDrive/Work/FACET2/GPT';% CE Change 6/2023
    gptloc='/usr/local/bin/';


    %NMM: Beam specs
    nmacro=1e5; % # of macro-particles. 1e4 is rough scan, 1e5 is point of diminishing returns, per CE
    Q0=[0.5e-9 1.6e-9]; % Initial bunch charge
    %laserR=2.75e-3; %0.5e-3; %2.75e-3; % cut radius, distribution is rms with 2XlaserR sigma (2.68) C3=1.88
    laserSigmaR = 1.7e-3; %NMM new
    laserCutR = 2.7e-3;
    dt=11e-12; %10.75e-12; % no LH=11e-12 LH,500=11.62e-12 % Offset between drive and witness beams
    sigt=[2.2e-12 2.2e-12]; % Laser pulse duration per bunch (FWHM in s) [2e-12 5.25e-12] % NMM 2023-09-07, was 3.84 ps
    gptTimeOffset = 5e-12; %NMM


    %NMM: Gun specs
    gunEz=120; % on-crest gun accelerating gradient (MV/m)
    %gunPhase=280; % gun RF phase offset (degrees)
    gunSolB=0.252; % Peak gun solenoid field strength (T) [Integrated: 1= 0.1936 T.m]
    nsigcut=3; % Cut threhold for gaussian random number generator
    gunFrequency = 2.856e9; %NMM, Hz
    %Matlab is trash and you apparently can't reference variables here...
    %moving these new definitions down
    %timeOffsetPhaseCorrection = (obj.gptTimeOffset*obj.gunFrequency*360); %NMM, degrees
    schottkyPhase = 30; %NMM; This is the phase as we typically refer to it, i.e. zero crossing at zero. GPT phase assumes maximum launch field at "zero"
    %gunPhase = 270 + obj.schottkyPhase + obj.timeOffsetPhaseCorrection; %NMM Translated to GPT phase
    gunPhase;

    %NMM: Linac specs
    L0a_phase=154; % L0a phase (degrees) [peak accelerating phase] %133 %NMM: This was getting set to 154 in the other file; was 10 here
    L0_phaseOffset=0;
    L0a_V=18.5; %CE Change6/2023 to get 66 MeV at L0a exit%17.1; % L0a accelerating gradient (MV/m) % provides 59 MV of acceleration





    %NMM More uninteresting stuff
    BeamOut % Storage for beam at exit of simulation
    SimData
    StoreVec=true; % Store particle vectors
    nscan=30; % number of scan points (nscan^2 for 2d scans)
    scanpoints
    verbose=1;
    batchmethod='linscan'; % linscan | mltrain | rgen

  end

  properties
    Beam
    BeamSingle
  end

  properties(SetAccess=private)
    nnet
  end

  properties(Access=private)
    tjitterval
  end

  properties(Constant)
    datafiles={'rfgundata.gdf' 'soldata.gdf' 'bunchData.gdf' 'bunchDataSingle.gdf' 'slwake.gdf' 'INJL0.saveline'};
  end
  
  methods
    %NMM: I don't think submitFunc() is in use, but it cannot be deleted.
    %I've lobotomized it though
    function data=submitFunc(obj,iseed)

    end





    function MakeCathodeBeam(obj)
      % sigT in m (rms), Len in s (FWHM)
      for ibunch=1:length(obj.sigt)
        %sigT=obj.laserR;
        Len=obj.sigt(ibunch);
        sigcut=obj.nsigcut;
        temitmult=1; % 1 for LCLS gun (0.9 um-rad /mm)
        npart=obj.nmacro;
        KE=1; % eV
        if ibunch==1
          obj.BeamSingle=CreateBlankBeam(1, length(obj.sigt), KE*1e-9, 1);
          obj.BeamSingle.Bunch.Q=obj.Q0;
        end
        
        % Set Gaussian width fraction of cut radius
        %transcut=sigT;
        %sigT=sigT*2;
        
        %xvals=randn(npart*10,1).*sigT;
        %yvals=randn(npart*10,1).*sigT;
        %figure; title('Initial');
        %histogram(xvals);

        %NMM update
        %For usual cases npart*10 is plenty, but for very very small
        %truncations going to 100 helps
        xvals = randn(npart*10,1).*obj.laserSigmaR;
        yvals = randn(npart*10,1).*obj.laserSigmaR;

        R=sqrt(xvals.^2+yvals.^2);
        xvals(R>obj.laserCutR)=[];
        %figure;
        %histogram(xvals); title('after transcut');
        xvals=xvals(1:npart);
        yvals(R>obj.laserCutR)=[];
        yvals=yvals(1:npart);
        %figure;
        %histogram(xvals); title('after sampling')
        gamma=1+(KE/1e9)/obj.emass;
        beta=sqrt(1-gamma^-2);
        v=obj.clight*beta;
        Len=Len/(2*sqrt(2*log(2)));
        zvals=randn(npart,1).*Len.*v;
        px=randnt(sigcut,npart,1).*0.0009*temitmult;
        py=randnt(sigcut,npart,1).*0.0009*temitmult;
        pz=abs(randnt(sigcut,npart,1).*2.9e-6)+0.001978;
%         xvals=xvals+randnt(1)*std(xvals)*obj.laserjitter(1); % do this with offsets in GPT.in file instead
%         yvals=yvals+randnt(1)*std(yvals)*obj.laserjitter(2);
        zoff=v*((ibunch-1)*obj.dt);
        if ibunch==1
          B=CreateBlankBeam(1, npart, KE*1e-9, 1);
          X=B.Bunch.x;
          X(1:5,:)=[xvals'; atan(px./pz)'; yvals'; atan(py./pz)'; zvals'];
          B.Bunch.Q=ones(1,npart).*(obj.Q0(ibunch)/npart);
        else
          X=[xvals'; atan(px./pz)'; yvals'; atan(py./pz)'; zoff+zvals'; B.Bunch.x(6,:)];
          B.Bunch.Q=[B.Bunch.Q ones(1,npart).*(obj.Q0(ibunch)/npart)];
          B.Bunch.stop=zeros(size(B.Bunch.Q));
        end

        if ibunch==1
          B.Bunch.x=X;
        else
          B.Bunch.x=[B.Bunch.x X];
        end
        obj.Beam=B;
      end
    end





    function [Bout,xv,yv,zv,ev,emitXvec]=readData(obj,iseed,cmd)
      % - read in GPT tracking results at screen location and store in
      % Lucretia Beam format
      nscreen=1;
      g=load_gdf(sprintf('result_%d.gdf',iseed));
      scr=flip(find(arrayfun(@(x) isfield(g(x).p,'position'),1:length(g))));
      if exist('cmd','var') && isequal(cmd,'single')
        BeamIn=obj.BeamSingle;
      else
        BeamIn=obj.Beam;
      end
      if ~isempty(scr)
        Bout=cell(1,nscreen);
        for iscr=1:nscreen
          id=(g(scr(iscr)).d.ID); % particle order in GPT file entries
          % - flag any missing missing particles and set these as stopped
          % particles in Lucretia beam definition data structure
          goodray=ismember(1:length(BeamIn.Bunch.Q),id) ;
          % - make cell array of output beams in Lucretia format
          if exist('cmd','var') && isequal(cmd,'single')
            Bout{iscr}=obj.BeamSingle;
          else
            Bout{iscr}=obj.Beam;
          end
          xang=atan(g(scr(iscr)).d.Bx./g(scr(iscr)).d.Bz);
          yang=atan(g(scr(iscr)).d.By./g(scr(iscr)).d.Bz);
          gamma=g(scr(iscr)).d.G;
          beta=sqrt(1-(1./gamma.^2));
          v=beta.*obj.clight;
          z=g(scr(iscr)).d.t.*v-g(scr(iscr)).d.z;
          E=obj.emass.*sqrt(gamma.^2-1);
          xv=g(scr(iscr)).d.x;
          yv=g(scr(iscr)).d.y;
          Bout{iscr}.Bunch.x(:,goodray)=[xv'; xang'; yv'; yang'; z'; E'];
          Bout{iscr}.Bunch.stop(~goodray)=1;
          Bout{iscr}.Bunch.x(5,goodray)=Bout{iscr}.Bunch.x(5,goodray)-mean(Bout{iscr}.Bunch.x(5,goodray));
          Bout{iscr}.Bunch.x(:,~goodray)=0;
          %fprintf("In the top loop");
        end
      else
        Bout=[];
      end
      % Read snapshot locations into vectors
      if nargout>1
        xv=nan(length(g),length(obj.Beam.Bunch.Q)); yv=xv; zv=xv; ev=xv;

        %NMM new
        betaXvec = xv; %Create initial shape of the array
        betaZvec = xv;
        emitXvec = nan(length(g),1);

        for ig=1:length(g)
          if any(ig==scr)
            continue;
          end
          id=(g(ig).d.ID); % particle order in GPT file entries
          xv(ig,id)=g(ig).d.x;
          yv(ig,id)=g(ig).d.y;
          zv(ig,id)=g(ig).d.z;
          beta=sqrt(g(ig).d.Bx.^2+g(ig).d.By.^2+g(ig).d.Bz.^2);
          gamma=1./sqrt(1-beta.^2);
          ev(ig,id)=gamma.*obj.emass;

          %2023-09-21 NMM new; sanity checks passed
    
          betaXvec(ig,id) = g(ig).d.Bx;
          betaZvec(ig,id) = g(ig).d.Bz;
          %fprintf("In the bottom loop");

          meanXSq = mean((g(ig).d.x).^2);
          meanXpSq = mean((g(ig).d.Bx./beta).^2);
          meanXXp = mean((g(ig).d.x).*(g(ig).d.Bx./beta));
          emitX = sqrt(meanXSq*meanXpSq - meanXXp^2);
          emitXvec(ig) = mean(beta)*mean(gamma)*emitX; %Vaguely sloppy, but I think it's fine
        end
      end
%       delete(sprintf('result_%d.gdf',iseed));
    end





    function writeBeamFile(obj,BeamIn,cmd)
      % - Write Lucrtia bunch particles out in GPT input file
      xGPT=BeamIn.Bunch.x;
      gamma=1+BeamIn.Bunch.x(6,:)./0.511e-3;
      beta=sqrt(1-gamma.^-2);
      beta_x=beta.*sin(xGPT(2,:)); beta_y=beta.*sin(xGPT(4,:)); beta_z=sqrt(1-(beta_x./beta).^2-(beta_y./beta).^2).*beta;
      data=struct; data.d.x=xGPT(1,:); data.d.y=xGPT(3,:); % GPT data structure for pos coordinates
      % Convert z distribution into particle release times
      data.d.z=randn(1,length(BeamIn.Bunch.Q)).*0;
      v=beta*obj.clight;
      data.d.t=xGPT(5,:)./v;
      %Glen's old version
      %data.d.t=data.d.t-min(data.d.t);

      %NMM new version with fixed offset
      data.d.t=data.d.t+obj.gptTimeOffset;

      data.d.Bx=beta_x'; data.d.By=beta_y'; data.d.Bz=beta_z'; data.d.nmacro=BeamIn.Bunch.Q./obj.eQ';
      data.d.G=1./sqrt(1-(beta_x.^2+beta_y.^2+beta_z.^2));
      if exist('cmd','var') && isequal(cmd,'single')
        save_struct_to_gdf_file(fullfile(obj.srcdir,'bunchDataSingle.gdf'), data); % save Matlab data file to GPT .gdf file
      else
        save_struct_to_gdf_file(fullfile(obj.srcdir,'bunchData.gdf'), data); % save Matlab data file to GPT .gdf file
      end
    end






    function writeRunFile(obj,iseed,cmd)

      timeOffsetPhaseCorrection = -1*(obj.gptTimeOffset*obj.gunFrequency*360); %NMM, degrees
      fprintf('Time offset phase correction [degrees]: %.2f\n', timeOffsetPhaseCorrection);
      obj.gunPhase = 270 + obj.schottkyPhase + timeOffsetPhaseCorrection; %NMM Translated to GPT phase
      fprintf('GPT phase [degrees]: %.2f\n', obj.gunPhase);

      % Generate GPT lattice and run file
      % - require beam to have been setup
      if isempty(obj.Beam)
        error('No beam defined');
      end

      fid=fopen(sprintf('GPT_%d.in',iseed),'w'); % open txt file for generating GPT run instructions

      fprintf(fid,'m=me ;\n'); % electron mass
      fprintf(fid,'q=qe ;\n'); % electron charge

      fprintf(fid,'setfile("beam","%s") ;\n',fullfile(pwd,'bunchData.gdf')) ; % the bunch data file
     
      
      fprintf(fid,'Spacecharge3Dmesh("Cathode");\n'); % Cathode

      fprintf(fid,'map1D_TM("wcs","z",0,"%s","z","Ez",%g,%g,1.79447772373049e10) ;\n',fullfile(obj.runDirectory,'rfgundata.gdf'),obj.gunEz*1e6/2.4523,deg2rad(obj.gunPhase)); % RF Gun % NMM: 1.79447772373049e10 is S-band in radians/sec
      fprintf(fid,'map1D_B("wcs","z",0,"%s","z","Bz",%g) ;\n',fullfile(obj.runDirectory,'soldata.gdf'),obj.gunSolB); % Gun solenoid





      %NMM: New if statement to stop before linacs
      if true
        gamma=0.006/0.511e-3; % Design for 6 MeV
        rfL=3.095244; %3.0429; 3.095244
        rfCenter = 2.58;% CE edit = 2.1 - Glen's rfCenter = 2.5798 was givinng a large kick to the beam
        a0=0.1; %NMM, per GPT manual, this apparently sets the "initial attenuation constant"

        P=(obj.L0a_V*1e6)^2/(2*50e6*a0); %NMM, some more magic numbers in here. Broadly looks to be converting a voltage, L0a_V, to a power

        endL0A=4.1274; 

        dpha=0; % NMM: No jitter
      
        fprintf(fid,'trwlinac("wcs","z",%g,%g,54e6,%g,%g,%g,0,%g,1.79447772373049e10,%g) ;\n',rfCenter,a0,P,P,gamma,deg2rad(obj.L0a_phase+dpha),rfL) ; % L0a 6/18/21 measured to give 59 MeV gain  % NMM: 1.79447772373049e10 is S-band in radians/sec, 54e6 is presumably the shunt impedance
      
        fprintf(fid,'wakefield("wcs","z",%g,%g,%g,"%s","z","","","Wz") ;\n',rfCenter,rfL,50,fullfile(obj.runDirectory,'slwake.gdf')); % Wakefield %NMM: The "50" evidently is a "fringe parameter" which, delightfully, does not appear to be defined in the GPT manual

        fprintf(fid,'screen("wcs","I",%g);\n',endL0A); % Beam output location
        fprintf(fid,'zminmax("wcs","z",0,-0.001,%.8g) ;\n',endL0A+0.001); % Cut backwards going particles

      else
        screenPosition = 1.2;
        fprintf(fid,'screen("wcs","I",%g);\n',screenPosition); % Beam output location
        fprintf(fid,'zminmax("wcs","z",0,-0.001,%.8g) ;\n',screenPosition+0.001); % Cut backwards going particles

      end
      

      if obj.StoreVec
        fprintf(fid,'snapshot(%g,%g,%g) ;\n',0,15e-9,5.0e-12);
      end

      fclose(fid);

    end

  end

end

