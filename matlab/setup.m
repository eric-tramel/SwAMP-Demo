function setup
    newline = sprintf('\n');
    src = regexprep(ls('src/*/*.c'), '[ \t\n]+', ' ');
    to_remove = regexprep(ls('src/*/c*.c'), '[ \t\n]+', ' ');
    src_excl = strjoin( setdiff(strsplit(src), strsplit(to_remove)) )

    cxx = 'mex';
    cflags = '-g -largeArrayDims -I/usr/local/include -L/usr/local/lib -ldl CFLAGS="\$CFLAGS -std=c99"';

    fprintf('Compiling SwAMP... ')
    eval([cxx ' ' cflags ' -output bin/swamp src/swamp.c ' src_excl]);
    % eval([cxx ' ' cflags ' -output bin/swgamp src/swgamp.c ' src_excl]);   
    % eval([cxx ' ' cflags ' -lgsl -lgslcblas -output bin/cswgamp src/cswgamp.c ' src]);   

    % We probably need to check if GSL is available before executing this...
    
    fprintf('Done!\n')
end
