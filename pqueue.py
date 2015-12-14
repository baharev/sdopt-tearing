# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
import itertools
from heapq import heappop, heappush, heapify

# Based on: 
# docs.python.org/2/library/heapq.html#priority-queue-implementation-notes

class PriorityQueue:
    
    def __init__(self):
        # entry: [priority, count from counter, task]
        self.heap = []                     # list of entries arranged in a heap
        self.entry_finder = {}             # mapping of tasks to entries
        self.counter = itertools.count()   # unique sequence count
    
    def __getitem__(self, task):
        return self.entry_finder[task][0]

#     def __iter__(self):
#         return iter(self.entry_finder)

    def populate(self, iterable):
        assert not self.entry_finder, 'The queue is supposed to be empty!'
        heap = [ ]
        for task, priority in iterable:
            entry = [priority, next(self.counter), task]
            heap.append(entry)
            self.entry_finder[task] = entry
        heapify(heap)
        self.heap = heap        

    def __len__(self):
        return len(self.entry_finder)

    def peekitem(self):
        while self.heap and self.heap[0][2] is None:
            heappop(self.heap)
        if self.heap:
            priority, _, task = self.heap[0]
            return (task, priority)
        raise KeyError('the priority queue is empty')
    
    def __setitem__(self, task, priority):
        'Add a new task or update the priority of an existing task.'
        if task in self.entry_finder:
            entry = self.entry_finder.pop(task)
            entry[-1] = None
        count = next(self.counter)
        entry = [priority, count, task]
        self.entry_finder[task] = entry
        heappush(self.heap, entry)
    
    def __delitem__(self, task):
        'Mark an existing task as removed (None). Raise KeyError if not found.'
        entry = self.entry_finder.pop(task)
        entry[-1] = None
    
    def popitem(self):
        'Remove and return the lowest priority task. Raise KeyError if empty.'
        while self.heap:
            priority, _, task = heappop(self.heap)
            if task is not None:
                del self.entry_finder[task]
                return priority, task
        raise KeyError('pop from an empty priority queue')
