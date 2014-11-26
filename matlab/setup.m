function setup
    newline = sprintf('\n');
    src = regexprep(ls('src/*/*.c'), '[ \t\n]+', ' ');
    to_remove = regexprep(ls('src/*/c*.c'), '[ \t\n]+', ' ');
    src_excl = strjoin( setdiff(strsplit(src), strsplit(to_remove)) );

    cxx = 'mex';
    cflags = '-O -largeArrayDims -I/usr/local/include -L/usr/local/lib -ldl CFLAGS="\$CFLAGS -std=c99"';

    fprintf('Compiling SwAMP... ')
    eval([cxx ' ' cflags ' -output bin/swamp src/swamp.c ' src_excl]);
    eval([cxx ' ' cflags ' -output bin/swgamp src/swgamp.c ' src_excl]);   
    eval([cxx ' ' cflags ' -lgsl -lgslcblas -output bin/cswgamp src/cswgamp.c ' src]);   

    % We probably need to check if GSL is available before executing this...
    

    % We should probably check and see if we have any lingering output
    % files, just in case. This seems to be an issue with building OSX.
    eval('delete *.o');

    fprintf('Done!\n')
end
