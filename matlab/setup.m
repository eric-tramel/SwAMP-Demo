function setup
    newline = sprintf('\n');
    src = regexprep(ls('src/*/*.c'), '[ \t\n]+', ' ');
    cxx = 'mex';
    cflags = '-largeArrayDims';

    fprintf('Compiling SwAMP... ')
    eval([cxx ' ' cflags ' -output bin/swamp src/swamp.c ' src])
    eval([cxx ' ' cflags ' -output bin/swgamp src/swgamp.c ' src]);
    fprintf('Done!\n')
end
