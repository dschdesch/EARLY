function Y = subsref(P, S);
% ActionButton/SUBSREF - SUBSREF for ActionButton objects
%
%   A(5) is the 5th element of ActionButton array A.
%   A('foo') is the element of array P whose name is Foo (exact match), [] if
%   no element is named Foo.
%   A('Foo', 'Goo', ..) is the same as [A('Foo') A('Goo') ..]
%
%   A('Foo*') returns the subarray of array A whose names start with Foo.
%
%   P.Foo is the the same as Foo(P) provided that Foo is a method for tsig
%   objects.
%
%   See also tsig, "methods tsig".

if length(S)>1, % use recursion from the left
    y = subsref(P,S(1));
    Y = subsref(y,S(2:end));
    return;
end
%----single-level subsref from here (i.e., S is scalar) ---
switch S.type,
    case '()', % array element
        if iscellstr(S.subs), % look for named object
            nam = S.subs;
            Y = [];
            if isequal(1,numel(nam)) && ~isempty(nam{1}) && isequal('*', nam{1}(end)), % wildcards
                imatch = strmatch(nam{1}(1:end-1),{P.Name});
                Y = P(imatch);
            else, % exact matches
                for ii=1:numel(nam),
                    imatch = strmatch(nam{ii},{P.Name},'exact');
                    if ~isempty(imatch), Y = [Y P(imatch)]; end;
                end
            end
            if isempty(Y),
                error(['Array of ' class(P) ' object does not contain elements with the name(s) as specified in parentheses.'])
            end
        else,
            Y = builtin('subsref',P,S);
        end
    case '.',
        Field = S.subs; % field name
        [Method, Mess] = keywordMatch(Field, methods(P), 'ActionButton field name or method');
        if isempty(Mess), 
            Y = feval(Method, P);
        else, % try fieldname
            [Fldname, Mess] = keywordMatch(Field, fieldnames(P), 'ActionButton field name or method');
            if isempty(Mess),
                Y = P.(Fldname);
            end
        end
        error(Mess);
    case '{}',
        error('ActionButton objects do not allow subscripted reference using braces {}.');
end % switch/case
