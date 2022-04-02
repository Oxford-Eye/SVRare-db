# given intervals, produce groups of intervals, where each member of a group overlaps with each other.
import numpy as np
import copy
from scipy.spatial import distance
from sklearn.cluster import DBSCAN


class Interval_base(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.size = end - start
        self.groups = []
        self.name = f'{chrom}-{start}-{end}'

    @staticmethod
    def overlap(this, other):
        if this.chrom == other.chrom:
            overlap_size = max(this.end, other.end) - \
                min(this.start, other.start)
            total_size = this.size + other.size
            if total_size > overlap_size:
                return True
        return False

    @staticmethod
    def get_distance(this, other):
        # 1 - (shared size / merged size)
        if this.chrom == other.chrom:
            distance = 1 - (min(this.end, other.end) - max(this.start, other.start)) / \
                (max(this.end, other.end) - min(this.start, other.start))

            return distance
        return None

    def __eq__(self, other):
        return ((self.chrom, self.start, self.end) == (other.chrom, other.start, other.end))

    def __lt__(self, other):
        if self.chrom < other.chrom:
            return True
        if self.chrom == other.chrom:
            if self.start < other.start:
                return True
            if self.start == other.start:
                if self.end < other.end:
                    return True
        return False

    def __le__(self, other):
        if self.chrom < other.chrom:
            return True
        if self.chrom == other.chrom:
            if self.start < other.start:
                return True
            if self.start == other.start:
                if self.end <= other.end:
                    return True
        return False

    def __gt__(self, other):
        if self.chrom > other.chrom:
            return True
        if self.chrom == other.chrom:
            if self.start > other.start:
                return True
            if self.start == other.start:
                if self.end > other.end:
                    return True
        return False

    def __ge__(self, other):
        if self.chrom > other.chrom:
            return True
        if self.chrom == other.chrom:
            if self.start > other.start:
                return True
            if self.start == other.start:
                if self.end >= other.end:
                    return True
        return False

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    @staticmethod
    def group(intervals):
        groups = []
        sorted_intervals = sorted(intervals)
        new_group = Group()
        for interval in sorted_intervals:
            check_interval = new_group.does_interval_belong_to_group(interval)
            if check_interval != False:
                new_group.add_interval(interval)
            else:
                groups.append(new_group)
                new_group = Group()
                new_group.add_interval(interval)
        groups.append(new_group)
        group_n = 0
        for interval in sorted_intervals:
            while group_n < len(groups) and interval in groups[group_n].intervals:
                group_n += 1
            while group_n < len(groups) and groups[group_n].does_interval_belong_to_group(interval):
                groups[group_n].add_interval(interval)
                group_n += 1
        return groups


class Group(object):
    def __init__(self):
        self.intervals = []
        self.chrom = None
        self.loose_start = None
        self.loose_end = None
        self.loose_name = None
        self.core_start = None
        self.core_end = None
        self.core_name = None
        self.core_interval = None
        self.loose_interval = None
        self.mean_start = None
        self.mean_end = None
        self.mean_interval = None
        self.mean_name = None
        self.stdev_start = None
        self.stdev_end = None

    def does_interval_belong_to_group(self, interval):
        # None: intervals empty, True: yes, False: no
        if not self.intervals:
            return None

        return Interval_base.overlap(interval, Interval_base(self.chrom, self.core_start, self.core_end))

    def add_interval(self, interval):
        if self.intervals and not self.does_interval_belong_to_group(interval):
            raise ValueError(
                f"interval: {interval} do not belong to group:{self} ")
        if self.intervals:
            self.loose_end = max(self.loose_end, interval.end)
            self.loose_start = min(self.loose_start, interval.start)
            self.core_start = max(
                self.core_start, interval.start)
            self.core_end = min(self.core_end, interval.end)
            self.loose_interval = Interval_base(
                self.chrom, self.loose_start, self.loose_end)
            self.core_interval = Interval_base(
                self.chrom, self.core_start, self.core_end)
            self.mean_start = int(round(np.mean(
                [i.start for i in self.intervals] + [interval.start])))
            self.stdev_start = np.std(
                [i.start for i in self.intervals] + [interval.start])
            self.mean_end = int(round(np.mean(
                [i.end for i in self.intervals] + [interval.end])))
            self.stdev_end = np.std(
                [i.end for i in self.intervals] + [interval.end])
            self.mean_interval = Interval_base(
                self.chrom, self.mean_start, self.mean_end)
        else:
            self.chrom = interval.chrom
            self.loose_end = self.core_end = self.mean_end = interval.end
            self.loose_start = self.core_start = self.mean_start = interval.start
            self.loose_interval = Interval_base(
                interval.chrom, interval.start, interval.end)
            self.core_interval = Interval_base(
                interval.chrom, interval.start, interval.end)
            self.mean_interval = Interval_base(
                interval.chrom, interval.start, interval.end)
            self.stdev_start = self.stdev_end = 0
        self.core_name = f'{self.chrom}-{self.core_start}-{self.core_end}'
        self.mean_name = f'{self.chrom}-{self.mean_start}-{self.mean_end}'
        self.intervals.append(interval)
        self.intervals.sort()

    def produce_distance_matrix(self, method=None):
        # use start and size to produce distance matrix
        if method is None:
            method = interval_distance_for_cdist
        A = [(i.start, i.size) for i in self.intervals]

        return distance.cdist(A, A, method)

    def mean_distance_to_core_interval(self):
        distances = [Interval_base.get_distance(self.core_interval, i)
                     for i in self.intervals]
        return {'mean': np.mean(distances), 'stdev': np.std(distances)}

    def mean_distance_to_mean_interval(self):
        distances = [Interval_base.get_distance(self.mean_interval, i)
                     for i in self.intervals]
        return {'mean': np.mean(distances), 'stdev': np.std(distances)}

    def cluster(self, distance):
        # making sub-groups with DBSCAN
        groups = {}
        D = self.produce_distance_matrix(
            method=interval_distance_for_cdist)
        C = DBSCAN(eps=distance, min_samples=2, metric='precomputed').fit(D)
        for group_index in set(C.labels_):
            if group_index == -1:
                for interval_index in [ind for ind in range(len(self.intervals)) if C.labels_[ind] == group_index]:
                    interval = self.intervals[interval_index]
                    group = Group()
                    group.add_interval(interval)
                    groups[group.core_name] = group
            else:
                intervals = [self.intervals[interval_index] for interval_index in range(
                    len(self.intervals)) if C.labels_[interval_index] == group_index]
                group = Group()
                for interval in intervals:
                    group.add_interval(interval)
                groups[group.core_name] = group

        return groups

    def normalise(self):
        pass

    def __eq__(self, other):
        return (self.core_interval == other.core_interval)

    def __lt__(self, other):
        return (self.core_interval < other.core_interval)

    def __le__(self, other):
        return (self.core_interval <= other.core_interval)

    def __gt__(self, other):
        return (self.core_interval > other.core_interval)

    def __ge__(self, other):
        return (self.core_interval >= other.core_interval)

    def __repr__(self):
        result = f'loose interval: {self.loose_interval}\n'
        result += f'core interval: {self.core_interval}\n'
        result += f'mean interval: {self.mean_interval}\n'
        result += f'stdev start: {self.stdev_start}\n'
        result += f'stdev end: {self.stdev_end}\n'
        result += f'intervals: {self.intervals}\n'
        result += f''
        return result


def interval_distance_for_cdist(u, v):
    # customised distance measure
    # 1 - (shared size / merged size)
    distance = 1 - (min(u[0]+u[1], v[0]+v[1]) - max(u[0], v[0])) / \
        (max(u[0]+u[1], v[0]+v[1]) - min(u[0], v[0]))
    return distance
