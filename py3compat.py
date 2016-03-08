# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>

# For some reason, it gives an annoying error in the IDE
from six.moves import zip_longest as izip_longest
from six.moves import zip as izip
from six.moves import xrange as irange
from six.moves import filter as ifilter
from six.moves import map as imap
from six.moves import cStringIO as StringIO
from six.moves import cPickle 
cPickle_loads = cPickle.loads # from string
cPickle_dumps = cPickle.dumps # to string
cPickle_load  = cPickle.load  # from file
cPickle_dump  = cPickle.dump  # to file
cPickle_HIGHEST_PROTOCOL = cPickle.HIGHEST_PROTOCOL
