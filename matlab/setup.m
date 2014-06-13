function setup
    newline = sprintf('\n');
    src = regexprep(ls('src/*/*.c'), '[ \t\n]+', ' ');
    cxx = 'mex';
    cflags = '-largeArrayDims';

    fprintf('Compiling SwAMP... ')
    eval([cxx ' ' cflags ' -o bin/swamp src/swamp.c ' src])
    eval([cxx ' ' cflags ' -o bin/swgamp src/swgamp.c ' src]);
    fprintf('Done!\n')
end
