function compile_test_debug

if ~strcmp(mexext,'mexmaci64')
  error('This tool has to be called from Mac OSX 64bit');
%  return
end

d0 = load('debug_mexmaci64.mat');

mexext_str = {'mexw64','mexa64','mexglx','mexw32'};

for i = 1:length(mexext_str)
  debugname = ['debug_' mexext_str{i} '.mat'];
  d = load(debugname);
  n = length(d0.d);
  fprintf('Absolute difference between maci64 and %s:\n',mexext_str{i});
  for j=1:n
    df = d0.d{j} - d.d{j};
    fprintf('%d: %f\n',j,max(df(:)));
  end
  j = j + 1;
  df = d0.CS.faces - d.CS.faces;
  fprintf('%d: faces %f\n',j,max(df(:)));
  df = d0.CS.vertices - d.CS.vertices;
  fprintf('%d: vertices %f\n',j,max(df(:)));
end