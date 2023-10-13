%addpath('/home/glen/OneDrive/Work/FACET2/GPT')% for save_struct_to_gdf_file.m and other functions
%% Track with GPT to the end of L0a
FI = GPT_FacetInjector_NMM();




%Write GPT beam and run files then run GPT
FI.MakeCathodeBeam;
FI.writeBeamFile(FI.Beam);
%FI.writeBeamFile(FI.BeamSingle,'single');
iseed = 1;
FI.writeRunFile(iseed);
sid=system(sprintf('gpt -o %s %s > /dev/null 2>&1',sprintf('result_%d.gdf',iseed),sprintf('GPT_%d.in',iseed)));


%% Read the simulation data
[Bout,xv,yv,zv,ev]=FI.readData(1);
beamImage(Bout{1})
[nx,ny,nz]=GetNEmitFromBeam(Bout{1},1);
[nx90,ny90,nz90,~,~]=GetNEmit90FromBeam(Bout{1})


figure; %NMM: This command prevents this new plot from being added as a subplot
plot(mean(zv,2,'omitnan'),1000*mean(ev,2,'omitnan')) %NMM: Plot average energy vs average z