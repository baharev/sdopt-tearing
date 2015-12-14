# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>

# For some reason, it gives an annoying error in the IDE
from six.moves import zip_longest as izip_longest
from six.moves import zip as izip
from six.moves import xrange as irange
from six.moves import cStringIO as StringIO
from six.moves import cPickle 
cPickle_loads = cPickle.loads
cPickle_dumps = cPickle.dumps 
cPickle_HIGHEST_PROTOCOL = cPickle.HIGHEST_PROTOCOL
