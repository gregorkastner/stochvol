def FlagsForFile (filename, **kwargs):
    return {
            'flags': ['-x', 'c++', '-std=c++11', '-I', '/usr/share/R/include',
                '-I', '/home/dhosszejni/include',
                '-I', '/home/dhosszejni/R/x86_64-pc-linux-gnu-library/3.4/Rcpp/include',
                '-I', '/home/dhosszejni/R/x86_64-pc-linux-gnu-library/3.4/RInside/include',
                '-I', '/home/dhosszejni/R/x86_64-pc-linux-gnu-library/3.4/RcppArmadillo/include',
                '-I', 'inst/include',
                '-fopenmp=libomp', '-fpic',
                '-Wall', '-Wformat', '-Werror=format-security', '-Wdate-time',
                '-g', '-c'],
            }
