function dmd_run_quick_function(app)

fun_name = app.QuickFunctionsDropDown.Value;
args = app.FunArgsTextArea.Value;
% arg_idx = strfind(fun_name,'()');
if isempty(args{:})
    feval(string(fun_name),app);
else
    eval([fun_name '(' 'app,' char(args) ')'])
end
