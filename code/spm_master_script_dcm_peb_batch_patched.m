%---------------------------------------------------------------------------------------
%                                                            
%      .d8888b.           888                      
%     d88P  Y88b          888                      
%     Y88b.               888                      
%      "Y888b.    .d88b.  888888 888  888 88888b.  
%         "Y88b. d8P  Y8b 888    888  888 888 "88b 
%           "888 88888888 888    888  888 888  888 
%     Y88b  d88P Y8b.     Y88b.  Y88b 888 888 d88P 
%      "Y8888P"   "Y8888   "Y888  "Y88888 88888P"  
%                                         888      
%                                         888      
%                                         888                                               
%                                                                        
%---------------------------------------------------------------------------------------
% Set up the MATLAB & SPM environment with necessary paths and variables

%---------------------------------------------------------------------------------------
% STEP 1 
%---------------------------------------------------------------------------------------

% Add SPM12 to MATLAB Path
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
rmpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
addpath('/home/py01/Projects/spm12') % Contains modified spm_dcm_fit with parfor=true

%---------------------------------------------------------------------------------------
% STEP 2 
%---------------------------------------------------------------------------------------

% Initialize SPM
spm('asciiwelcome');
spm_jobman('initcfg'); % Allows batch operations inside a script
spm('defaults','EEG');
spm_get_defaults('cmdline',true);

%---------------------------------------------------------------------------------------
% STEP 3 
%---------------------------------------------------------------------------------------

% Specify root working directory 
base_dir = '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/';
addpath(fullfile(base_dir, 'code')) % Add scripts & functions to workspace

%---------------------------------------------------------------------------------------
%                                                 
%     8888888b.   .d8888b.  888b     d888 
%     888  "Y88b d88P  Y88b 8888b   d8888 
%     888    888 888    888 88888b.d88888 
%     888    888 888        888Y88888P888 
%     888    888 888        888 Y888P 888 
%     888    888 888    888 888  Y8P  888 
%     888  .d88P Y88b  d88P 888   "   888 
%     8888888P"   "Y8888P"  888       888
%
%---------------------------------------------------------------------------------------
% Configure DCM model with data and options

% Name of this DCM model (Used for each subject's fits, and downstream by GCM/PEB)
name = 'DCM_MEG_ECD_OFA-FFA_Full';
DCM.name = name;

%---------------------------------------------------------------------------------------
% STEP 1: Setup analysis options
%---------------------------------------------------------------------------------------

% Specify modality
DCM.xY.modality = 'MEG';

% Set up DCM analysis type
DCM.options.analysis = 'ERP';   % Analyze evoked responses
DCM.options.model    = 'ERP';   % Neuronal temporal model: Extended Jansen-Rit model
DCM.options.spatial  = 'ECD';   % Spatial observation model: ECD

% Set up preprocessing parameters and analysis options
DCM.options.Nmodes   = 8;       % Number of modes of Leadfield for data selection
DCM.options.CVA      = 0;       % Optimize modes of Leadfield
DCM.options.h        = 1;       % Number of DCT components for detrending (1 is mean)
DCM.options.han      = 1;       % Hanning Window Taper
DCM.options.onset    = 64;      % Selection of onset (prior mean) for input stimulus
DCM.options.dur      = 16;      % Duration of onset (prior sd) for input stimulus
DCM.options.D        = 1;       % Downsampling (decimation of time series by a factor)
DCM.options.multiC   = 0;       % Multiple input vectors for multiple stimuli
DCM.options.location = 0;       % Optimize dipole locations
DCM.options.symmetry = 1;       % Lock orientation of dipoles across hemispheres

%---------------------------------------------------------------------------------------
% STEP 2: Setup data & design
%---------------------------------------------------------------------------------------

% Specify data of interest
DCM.options.trials   = [1 2 3]; % Index of ERPs within ERP/ERF file
DCM.options.Tdcm     = [0 500]; % Peri-stimulus time to be modelled

% Specify between-condition trial effects
contrasts = [1 1 0]'; % Face Perception: Faces (Famous + Unfamiliar) - Scrambled
DCM.xU.X = contrasts; % Orientation is N_trials x N_contrasts
DCM.xU.name = {'Face Perception'};

%--------------------------------------------------------------------------
% STEP 3: Setup observation model
%--------------------------------------------------------------------------

% Location priors for dipoles
locs  = {
    [-38, -86, -14], 'lOFA';
    [+36, -86, -10], 'rOFA';
        
    [-42, -56, -20], 'lFFA';
    [+42, -52, -14], 'rFFA';   
};

DCM.Lpos  = cat(1, locs{:,1})';
DCM.Sname = locs(:,2)';
Nareas    = length(locs);

%--------------------------------------------------------------------------
% STEP 4: Setup neuronal model
%--------------------------------------------------------------------------

% A Matrix: Forward connections
DCM.A{1} = [
%    lOFA rOFA lFFA rFFA
    [  0    0    0    0  ];   % lOFA
    [  0    0    0    0  ];   % rOFA
    [  1    0    0    0  ];   % lFFA
    [  0    1    0    0  ];   % rFFA    
];

% A Matrix: Backward connections
DCM.A{2} = [
%    lOFA rOFA lFFA rFFA
    [  0    0    1    0  ];   % lOFA
    [  0    0    0    1  ];   % rOFA
    [  0    0    0    0  ];   % lFFA
    [  0    0    0    0  ];   % rFFA
];

% A Matrix: Lateral connections
DCM.A{3} = [
%    lOFA rOFA lFFA rFFA
    [  0    1    0    0  ];   % lOFA
    [  1    0    0    0  ];   % rOFA
    [  0    0    0    1  ];   % lFFA
    [  0    0    1    0  ];   % rFFA
];

% B Matrix: Modulation of connections
self_connections = eye(Nareas);
DCM.B{1} = double(DCM.A{1} | DCM.A{2} | DCM.A{3} | self_connections);

% C Matrix: Driving inputs
DCM.C = [1 1 0 0]';

% Save full model
dcm_full_file = fullfile(base_dir, 'fits', 'batch', strcat(DCM.name, '_Specification.mat'));
save(dcm_full_file, 'DCM')
DCM_Full = DCM;

%--------------------------------------------------------------------------
% STEP 5: Specify reduced models, if any
%--------------------------------------------------------------------------

% Reduced model with modulation of only self connections
DCM.name = 'DCM_MEG_ECD_OFA-FFA_Self';
DCM.B{1} = self_connections;

% Save reduced model
dcm_self_file = fullfile(base_dir, 'fits', 'batch', strcat(DCM.name, '_Specification.mat'));
save(dcm_self_file, 'DCM')

%---------------------------------------------------------------------------------------
%                                                                                         
%     8888888888         888    d8b                        888            
%     888                888    Y8P                        888            
%     888                888                               888            
%     8888888   .d8888b  888888 888 88888b.d88b.   8888b.  888888 .d88b.  
%     888       88K      888    888 888 "888 "88b     "88b 888   d8P  Y8b 
%     888       "Y8888b. 888    888 888  888  888 .d888888 888   88888888 
%     888            X88 Y88b.  888 888  888  888 888  888 Y88b. Y8b.     
%     8888888888 88888P'  "Y888 888 888  888  888 "Y888888  "Y888 "Y8888
%
%---------------------------------------------------------------------------------------
% Estimate specified DCMs for all subjects using the batch interface

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Cell array of paths to DCM files
dcm_files = {dcm_full_file; dcm_self_file};

% String identifier for estimated GCM & PEB files
name_tag = 'MEG_ECD_OFA-FFA_Fit';

% Populate list of processed files as a column-order cell (N-files × 1)
% These files should contain forward models (with or without gain matrices)
files = dir(fullfile(base_dir, 'data', '**', 'maM*.mat'));
input_files = arrayfun(@(x) fullfile(x.folder, x.name), files, 'UniformOutput', false);

% Specify output directory
out_dir = {fullfile(base_dir, 'fits', 'batch')};

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(base_dir, 'code', 'batch_dcm_peb_estimation_job.m')};

% Initialize inputs to jobfile as a columnar cell array
% There are 6 inputs needed by this jobfile
inputs = cell(6, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1, 1} = input_files; % DCM for M/EEG: M/EEG datasets - cfg_files (Should be cell array, row-order)
inputs{2, 1} = dcm_files; % DCM for M/EEG: DCM files - cfg_files (Should be cell array, not string/char)
inputs{3, 1} = out_dir; % DCM for M/EEG: Directory - cfg_files (Should be cell array, not string/char)
inputs{4, 1} = out_dir; % DCM estimation: Directory - cfg_files (Should be cell array, not string/char)
inputs{5, 1} = name_tag; % DCM estimation: Name - cfg_entry ('GCM_' suffix will be appended)
inputs{6, 1} = name_tag; % Specify / Estimate PEB: Name - cfg_entry ('PEB_' suffix will be appended)

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job using parallel compute backend
%---------------------------------------------------------------------------------------

% Initialize Parallel Compute Pool (Example Instructions for CBU Cluster)
delete(gcp('nocreate')) % Shut down any existing pool
n_workers = length(input_files);
P=cbupool(n_workers, '--mem-per-cpu=4G --time=12:00:00 --nodelist=node-j10');
parpool(P, P.NumWorkers);

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
%
%      .d8888b.                                    888      
%     d88P  Y88b                                   888      
%     Y88b.                                        888      
%      "Y888b.    .d88b.   8888b.  888d888 .d8888b 88888b.  
%         "Y88b. d8P  Y8b     "88b 888P"  d88P"    888 "88b 
%           "888 88888888 .d888888 888    888      888  888 
%     Y88b  d88P Y8b.     888  888 888    Y88b.    888  888 
%      "Y8888P"   "Y8888  "Y888888 888     "Y8888P 888  888
%
%---------------------------------------------------------------------------------------
% Perform greedy search over full model space (B-matrix) with Bayesian model reduction

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Paths to estimated PEB and GCM structures
PEB_file = cellstr(fullfile(out_dir, sprintf('PEB_%s.mat', name_tag)));
GCM_file = cellstr(fullfile(out_dir, sprintf('GCM_%s.mat', name_tag)));

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(base_dir, 'code', 'batch_dcm_peb_search_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 2 inputs needed by this jobfile
inputs = cell(2, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1} = PEB_file;
inputs{2} = GCM_file;

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
%
%      .d8888b.                                                           
%     d88P  Y88b                                                          
%     888    888                                                          
%     888         .d88b.  88888b.d88b.  88888b.   8888b.  888d888 .d88b.  
%     888        d88""88b 888 "888 "88b 888 "88b     "88b 888P"  d8P  Y8b 
%     888    888 888  888 888  888  888 888  888 .d888888 888    88888888 
%     Y88b  d88P Y88..88P 888  888  888 888 d88P 888  888 888    Y8b.     
%      "Y8888P"   "Y88P"  888  888  888 88888P"  "Y888888 888     "Y8888  
%                                       888                               
%                                       888                               
%                                       888
%
%---------------------------------------------------------------------------------------
% Perform Bayesian Model Comparison (BMC) for Full vs reduced Self-only models

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Paths to estimated PEB and GCM structures
PEB_file = cellstr(fullfile(out_dir, sprintf('PEB_%s.mat', name_tag)));
GCM_file = cellstr(fullfile(out_dir, sprintf('GCM_%s.mat', name_tag)));

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(base_dir, 'code', 'batch_dcm_peb_bmc_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 2 inputs needed by this jobfile
inputs = cell(2, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1} = PEB_file;
inputs{2} = GCM_file;

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
%
%     888b     d888                               888888b.   888b     d888  .d8888b.  
%     8888b   d8888                               888  "88b  8888b   d8888 d88P  Y88b 
%     88888b.d88888                               888  .88P  88888b.d88888 888    888 
%     888Y88888P888  .d88b.  888d888 .d88b.       8888888K.  888Y88888P888 888        
%     888 Y888P 888 d88""88b 888P"  d8P  Y8b      888  "Y88b 888 Y888P 888 888        
%     888  Y8P  888 888  888 888    88888888      888    888 888  Y8P  888 888    888 
%     888   "   888 Y88..88P 888    Y8b.          888   d88P 888   "   888 Y88b  d88P 
%     888       888  "Y88P"  888     "Y8888       8888888P"  888       888  "Y8888P"
%
%---------------------------------------------------------------------------------------
% Since we observed significant modulation of between-region connections, we can zoom
% in and test for modulation of selective groups of connections like forward, backward
% and lateral connections. We demo this sequential hypothesis testing here using BMC.

%---------------------------------------------------------------------------------------
%
%     8888888888                                                  888 
%     888                                                         888 
%     888                                                         888 
%     8888888  .d88b.  888d888 888  888  888  8888b.  888d888 .d88888 
%     888     d88""88b 888P"   888  888  888     "88b 888P"  d88" 888 
%     888     888  888 888     888  888  888 .d888888 888    888  888 
%     888     Y88..88P 888     Y88b 888 d88P 888  888 888    Y88b 888 
%     888      "Y88P"  888      "Y8888888P"  "Y888888 888     "Y88888
%
%---------------------------------------------------------------------------------------
% Test for Modulation of Forward Connections

%---------------------------------------------------------------------------------------
% STEP 1: Define model space
%---------------------------------------------------------------------------------------

% Get full DCM specification
DCM = DCM_Full;

% Switch off Forward connections in B-matrix
DCM.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  1    1    1    0  ];   % lOFA
    [  1    1    0    1  ];   % rOFA
    [  0    0    1    1  ];   % lFFA
    [  0    0    1    1  ];   % rFFA
];

% Save GCM
dcm_forwardOff_file = fullfile(out_dir, 'DCM_ForwardOFF_Specification.mat');
save(dcm_forwardOff_file, 'DCM')

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Paths to DCM models and estimated PEB
PEB_file = cellstr(fullfile(out_dir, sprintf('PEB_%s.mat', name_tag)));
DCM_files = {dcm_full_file; dcm_forwardOff_file};

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(base_dir, 'code', 'batch_dcm_peb_bmc_modelspace_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 4 inputs needed by this jobfile
inputs = cell(4, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1} = input_files;
inputs{2} = DCM_files;
inputs{3} = out_dir;
inputs{4} = PEB_file;

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

% No forward connections from OFA to FFA are not modulated by Faces

%---------------------------------------------------------------------------------------
%
%     888888b.                     888                                         888 
%     888  "88b                    888                                         888 
%     888  .88P                    888                                         888 
%     8888888K.   8888b.   .d8888b 888  888 888  888  888  8888b.  888d888 .d88888 
%     888  "Y88b     "88b d88P"    888 .88P 888  888  888     "88b 888P"  d88" 888 
%     888    888 .d888888 888      888888K  888  888  888 .d888888 888    888  888 
%     888   d88P 888  888 Y88b.    888 "88b Y88b 888 d88P 888  888 888    Y88b 888 
%     8888888P"  "Y888888  "Y8888P 888  888  "Y8888888P"  "Y888888 888     "Y88888
%
%---------------------------------------------------------------------------------------
% Test for Modulation of Backward Connections

%---------------------------------------------------------------------------------------
% STEP 1: Define model space
%---------------------------------------------------------------------------------------

% Get full DCM specification
DCM = DCM_Full;

% Switch off Backward connections in B-matrix
DCM.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  1    1    0    0  ];   % lOFA
    [  1    1    0    0  ];   % rOFA
    [  1    0    1    1  ];   % lFFA
    [  0    1    1    1  ];   % rFFA
];

% Save DCM
dcm_backwardOff_file = fullfile(out_dir, 'DCM_BackwardOFF_Specification.mat');
save(dcm_backwardOff_file, 'DCM')

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Paths to DCM models and estimated PEB
PEB_file = cellstr(fullfile(out_dir, sprintf('PEB_%s.mat', name_tag)));
DCM_files = {dcm_full_file; dcm_backwardOff_file};

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(base_dir, 'code', 'batch_dcm_peb_bmc_modelspace_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 4 inputs needed by this jobfile
inputs = cell(4, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1} = input_files;
inputs{2} = DCM_files;
inputs{3} = out_dir;
inputs{4} = PEB_file;

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

% Both backward connections from FFA to OFA are modulated by Faces

%---------------------------------------------------------------------------------------
%
%     888               888                            888 
%     888               888                            888 
%     888               888                            888 
%     888       8888b.  888888 .d88b.  888d888 8888b.  888 
%     888          "88b 888   d8P  Y8b 888P"      "88b 888 
%     888      .d888888 888   88888888 888    .d888888 888 
%     888      888  888 Y88b. Y8b.     888    888  888 888 
%     88888888 "Y888888  "Y888 "Y8888  888    "Y888888 888
%
%---------------------------------------------------------------------------------------
% Test for Modulation of Lateral Connections

%---------------------------------------------------------------------------------------
% STEP 1: Define model space
%---------------------------------------------------------------------------------------

% Get full DCM specification
DCM = DCM_Full;

% Switch off Lateral connections in B-matrix
DCM.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  1    0    1    0  ];   % lOFA
    [  0    1    0    1  ];   % rOFA
    [  1    0    1    0  ];   % lFFA
    [  0    1    0    1  ];   % rFFA
];

% Save DCM
dcm_lateralOff_file = fullfile(out_dir, 'DCM_LateralOFF_Specification.mat');
save(dcm_lateralOff_file, 'DCM')

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Paths to DCM models and estimated PEB
PEB_file = cellstr(fullfile(out_dir, sprintf('PEB_%s.mat', name_tag)));
DCM_files = {dcm_full_file; dcm_lateralOff_file};

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(base_dir, 'code', 'batch_dcm_peb_bmc_modelspace_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 4 inputs needed by this jobfile
inputs = cell(4, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1} = input_files;
inputs{2} = DCM_files;
inputs{3} = out_dir;
inputs{4} = PEB_file;

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

% Bidirectional lateral connections between both OFA and FFA are modulated by Faces
% Follow-up: Test for OFA and FFA separately?

%---------------------------------------------------------------------------------------
%
%      .d8888b.           888  .d888 
%     d88P  Y88b          888 d88P"  
%     Y88b.               888 888    
%      "Y888b.    .d88b.  888 888888 
%         "Y88b. d8P  Y8b 888 888    
%           "888 88888888 888 888    
%     Y88b  d88P Y8b.     888 888    
%      "Y8888P"   "Y8888  888 888
%
%---------------------------------------------------------------------------------------
% After testing which between-region connections are modulated, we test which self
% connections are modulated by faces.

%---------------------------------------------------------------------------------------
% STEP 1: Define model space
%---------------------------------------------------------------------------------------

% Get full DCM specification
DCM = DCM_Full;

% Switch off Self connections in B-matrix
DCM.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    1    0  ];   % lOFA
    [  1    0    0    1  ];   % rOFA
    [  1    0    0    1  ];   % lFFA
    [  0    1    1    0  ];   % rFFA
];

% Save DCM
dcm_selfOff_file = fullfile(out_dir, 'DCM_SelfOFF_Specification.mat');
save(dcm_selfOff_file, 'DCM')

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Paths to DCM models and estimated PEB
PEB_file = cellstr(fullfile(out_dir, sprintf('PEB_%s.mat', name_tag)));
DCM_files = {dcm_full_file; dcm_selfOff_file};

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(base_dir, 'code', 'batch_dcm_peb_bmc_modelspace_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 4 inputs needed by this jobfile
inputs = cell(4, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1} = input_files;
inputs{2} = DCM_files;
inputs{3} = out_dir;
inputs{4} = PEB_file;

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

% No self connections in OFA and FFA are modulated by Faces

%---------------------------------------------------------------------------------------
%
%     8888888888                     d8b 888 d8b                   
%     888                            Y8P 888 Y8P                   
%     888                                888                       
%     8888888  8888b.  88888b.d88b.  888 888 888  .d88b.  .d8888b  
%     888         "88b 888 "888 "88b 888 888 888 d8P  Y8b 88K      
%     888     .d888888 888  888  888 888 888 888 88888888 "Y8888b. 
%     888     888  888 888  888  888 888 888 888 Y8b.          X88 
%     888     "Y888888 888  888  888 888 888 888  "Y8888   88888P'
%
%---------------------------------------------------------------------------------------
% Instead of testing for modulation of all self-connections versus none, we could ask
% whether any self-connection is being modulated at all. This involves multiple models
% with different combinations of self-connections being modulated. The model space
% therefore expands from just 2 models, like the previous BMC, to a much larger number
% of models. As an example here, we consider bilateral pairs of OFA and FFA, and test
% whether at least one self-connection was modulated by faces. These leads to a model
% space consisting of these four models:
%       1. Modulation of both OFA and FFA self-connections
%       2. Modulation of only OFA self-connections
%       3. Modulation of only FFA self-connections
%       4. Modulation of no self-connections
% The first three models encapsulate the hypothesis: are any self-connections modulated
% by faces? While the fourth model corresponds to the alternate hypothesis that no
% self-connections are modulated. Accordingly, we group these models into two families,
% and perform inference at the level of model families, rather than at the level of the
% individual models.

%---------------------------------------------------------------------------------------
% STEP 1: Setup model space and families
%---------------------------------------------------------------------------------------

% >> Family 1: Atleast one self-connection
% > Model 1 is the full model with both OFA and FFA self-connections
if isfield(DCM_Full, 'M')
    DCM_Full = rmfield(DCM_Full, 'M');
end
DCM_m1 = DCM_Full;

% > Model 2 has OFA self-connections only
% Get full DCM specification
DCM_m2 = DCM_m1;

% Switch off FFA self-connections in B-matrix
DCM_m2.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  1    1    1    0  ];   % lOFA
    [  1    1    0    1  ];   % rOFA
    [  1    0    0    1  ];   % lFFA
    [  0    1    1    0  ];   % rFFA
];

% > Model 3 has FFA self-connections only
% Get full DCM specification
DCM_m3 = DCM_m1;

% Switch off OFA self-connections in B-matrix
DCM_m3.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    1    0  ];   % lOFA
    [  1    0    0    1  ];   % rOFA
    [  1    0    1    1  ];   % lFFA
    [  0    1    1    1  ];   % rFFA
];

% >> Family 2: No self-connections
% > Model 4 has no self-connections
% Get full DCM specification
DCM_m4 = DCM_m1;

% Switch off all self-connections in B-matrix
DCM_m4.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    1    0  ];   % lOFA
    [  1    0    0    1  ];   % rOFA
    [  1    0    0    1  ];   % lFFA
    [  0    1    1    0  ];   % rFFA
];

% Define family-wise model space
models = {DCM_m1, DCM_m2, DCM_m3, DCM_m4};
families = [1, 1, 1, 2];

%---------------------------------------------------------------------------------------
% STEP 2: Perform BMR of model space nested under PEB
%---------------------------------------------------------------------------------------

% Bayesian Model Reduction (BMR) and comparison of models
[BMA, BMR] = spm_dcm_peb_bmc(PEB, models);
% Best model appears to be the one with modulation of only FFA self-connections
% This contradicts the previous BMC, where no self connections were modulated by Faces

% We now group models under families & consider each family as equally likely. This lets
% us pool model evidence within each family, and instead of picking a winning model, we
% perform inference about characteristics that define models grouped under that family.
% The characteristic here being modulation of self-connections, so models with any self
% connections i.e. models 1, 2 and 3, belong to the same family.

%---------------------------------------------------------------------------------------
% STEP 3: Compare families of models
%---------------------------------------------------------------------------------------

% Bayesian Model Selection (BMS) and averaging over families of models
[BMAf, fam] = spm_dcm_peb_bmc_fam(BMA, BMR, families, 'NONE');

% The family of models with modulation of at least one self-connection has marginally
% higher posterior probability than the family with no modulation of self-connections.
% There is not enough evidence in favor of modulation of self-connections by faces.

%---------------------------------------------------------------------------------------
%
%      .d8888b.                                     d8b          888                     
%     d88P  Y88b                                    Y8P          888                     
%     888    888                                                 888                     
%     888         .d88b.  888  888  8888b.  888d888 888  8888b.  888888 .d88b.  .d8888b  
%     888        d88""88b 888  888     "88b 888P"   888     "88b 888   d8P  Y8b 88K      
%     888    888 888  888 Y88  88P .d888888 888     888 .d888888 888   88888888 "Y8888b. 
%     Y88b  d88P Y88..88P  Y8bd8P  888  888 888     888 888  888 Y88b. Y8b.          X88 
%      "Y8888P"   "Y88P"    Y88P   "Y888888 888     888 "Y888888  "Y888 "Y8888   88888P'
%
%---------------------------------------------------------------------------------------
% We demonstrate the inclusion of covariates for 2nd-level (group) inference with PEB.
% The dataset includes ages of participants, and while we do not anticipate any effect
% of age on modulation of connections due to faces, we illustrate specification of age
% as a covariate in the PEB design matrix for inference.

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Define covariates, and assign appropriate labels
PEB_name = 'Age';
covariate_name = 'Age';
covariate_values = [31, 25, 30, 26, 23, 26, 31, 26, 29, 23, 24, 24, 25, 24, 30, 25]';
%covariate_name = 'Sex';
%covariate_values = [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0]'; % 1=Female

% Mean-center the covariate (Optional)
%covariate_values = covariate_values - mean(covariate_values);

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(base_dir, 'code', 'batch_dcm_peb_covariate_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 4 inputs needed by this jobfile
inputs = cell(4, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1, 1} = PEB_name; % Specify / Estimate PEB: Name - cfg_entry
inputs{2, 1} = GCM_file; % Specify / Estimate PEB: DCMs - cfg_files
inputs{3, 1} = covariate_name; % Specify / Estimate PEB: Name - cfg_entry
inputs{4, 1} = covariate_values; % Specify / Estimate PEB: Value - cfg_entry

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

spm_jobman('run', jobfile, inputs{:});

% Age has no effect on modulation of connections due to faces except for a small effect
% on rFFA self-connection.
