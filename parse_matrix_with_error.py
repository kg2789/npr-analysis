from uncertainties import ufloat_fromstr
import fileinput

# To be used for parsing the output of Greg's program and converting to a form
# which is easier to plug into Eigen

def parse_float_with_error(s):
    if s == "~0" or s == "0":
        return (0.0, 1e-10)
    f = ufloat_fromstr(s)
    value = f.nominal_value
    error = f.std_dev

    return (value, error)

values = []
errors = []

for line in fileinput.input():
    for num in line.rstrip().split():
        value, error = parse_float_with_error(num)
        values.append(value)
        errors.append(error)



print("Values:")
for val in values:
    print(val, end = ', ')

print()

print("Errors:")
for err in errors:
    print(err, end = ', ')

print()


