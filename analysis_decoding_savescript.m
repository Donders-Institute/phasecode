vararg = [];
% make filename
filename = [projectdir 'results/collapsephasebin/'];

switch contrast
  case 'congruent'
    filename = [filename, 'congruent/'];
  case 'attended'
    filename = [filename, 'attended/'];
  case 'unattended'
    filename = [filename, 'unattended/'];
end
filename = [filename, sprintf('sub%02d/', subj)];

if doparc
  filename = [filename, 'parc/'];
elseif strcmp(chansel_orig, 'eye')
  filename = [filename, 'eye/'];
elseif 0% exist('primal_P','var') && do_randphasebin
  filename = [filename, 'primalp/'];
end

if ~isdir(filename)
  mkdir(filename)
end

filename = [filename, sprintf('sub%02d_decoding', subj)];

if do_randphasebin
  filename = [filename, '_rand'];
end
if strcmp(contrast, 'attended') || strcmp(contrast, 'unattended')
  filename = [filename, sprintf('_hemi%d', hemi)];
end

if nbins>1
  filename = [filename, sprintf('_f%d', f)];
end

if doparc
  filename = [filename, sprintf('_%d', whichparc)];
end

if exist('randnr', 'var') && ~isempty(randnr)
  filename = [filename, sprintf('_%d', randnr)];
end

% save variables
settings = struct('correcttrials', do_correcttrials, 'contrast', contrast,'prewhiten', do_prewhiten>0, 'var', vararg, 'nperm',nperm, 'nrandperm', nrandperm);
if ~exist('primal', 'var') || do_randphasebin, primal=[];P=[]; end
rnumb=rng;
if dosave
  if isempty(tmpfilename)
    save(filename, 'accuracy','settings', 'primal_P', 'rnumb')
  else
    save(tmpfilename, 'accuracy', 'rnumb')
  end
end