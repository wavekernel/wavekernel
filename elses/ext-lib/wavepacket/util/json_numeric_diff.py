import json
import sys
import types

diff_float_tol = 0.01
filename1 = ''
filename2 = ''

def is_string_like(t):
    return t == types.StringType or t == types.UnicodeType

def path_to_str(path):
    res = ''
    for p in path:
        if is_string_like(type(p)):
            res += '["%s"]' % p
        elif type(p) == types.IntType:
            res += '[%d]' % p
        else:
            assert(False)
    return res

def diff(path, j1, j2):
    t1 = type(j1)
    t2 = type(j2)
    if t1 != t2:
        print ("%s:" % path_to_str(path)), t1, t2
    elif is_string_like(t1):
        diff_string(path, j1, j2)
    elif t1 == types.IntType:
        diff_int(path, j1, j2)
    elif t1 == types.FloatType:
        diff_float(path, j1, j2)
    elif t1 == types.BooleanType:
        diff_boolean(path, j1, j2)
    elif t1 == types.DictType:
        diff_dict(path, j1, j2)
    elif t1 == types.ListType:
        diff_list(path, j1, j2)
    elif t1 == types.NoneType:
        diff_none()
    else:
        print '%s: unknown type' % path_to_str(path), t1
    return

def diff_string(path, s1, s2):
    if s1 != s2:
        print '%s: "%s" "%s"' % (path_to_str(path), s1, s2)

def diff_int(path, i1, i2):
    if i1 != i2:
        print '%s: %d %d' % (path_to_str(path), i1, i2)

def diff_float(path, f1, f2):
    if abs(f1 - f2) > diff_float_tol:
        print '%s: %e %e' % (path_to_str(path), f1, f2)

def diff_boolean(path, b1, b2):
    if b1 != b2:
        print '%s: %r %r' % (path_to_str(path), b1, b2)

def make_keys_str(keys):
    keys_str = ''
    for k in keys:
        if keys_str == '':
            comma = ''
        else:
            comma = ', '
        keys_str += '%s"%s"' % (comma, k)
    return keys_str

def diff_dict(path, d1, d2):
    keys1 = set(d1.keys())
    keys2 = set(d2.keys())
    keys_only_in_1 = keys1.difference(keys2)
    keys_only_in_2 = keys2.difference(keys1)
    keys_in_both = keys1.intersection(keys2)
    if keys_only_in_1:
        print ('%s[%s]:' % (path_to_str(path), make_keys_str(keys_only_in_1))), ('<only in %s>' % filename1)
    if keys_only_in_2:
        print ('%s[%s]:' % (path_to_str(path), make_keys_str(keys_only_in_2))), ('<only in %s>' % filename2)
    for k in keys_in_both:
        diff(path + [k], d1[k], d2[k])

def diff_list(path, l1, l2):
    len1 = len(l1)
    len2 = len(l2)
    len_min = min(len1, len2)
    for i in range(len_min):
        diff(path + [i], l1[i], l2[i])
    if len1 > len2:
        if len1 == len2 + 1:
            str_index = "%d" % len2
        else:
            str_index = "%d - %d" % (len_min, len1 - 1)
        print ("%s[%s]:" % (path_to_str(path), str_index)), ('<only in %s>' % filename1)
    elif len2 > len1:
        if len2 == len1 + 1:
            str_index = "%d" % len1
        else:
            str_index = "%d - %d" % (len_min, len2 - 1)
        print ("%s[%s]:" % (path_to_str(path), str_index)), ('<only in %s>' % filename2)
    else:
        return

def diff_none():
    return

if __name__ == '__main__':
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print 'Usage: python %s <JSON> <JSON> [<float diff tol>]' % sys.argv[0]
        sys.exit(0)

    if len(sys.argv) == 4:
        diff_float_tol = float(sys.argv[3])

    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    with open(sys.argv[1]) as fp:
        json1 = json.load(fp)
    with open(sys.argv[2]) as fp:
        json2 = json.load(fp)

    diff([], json1, json2)
