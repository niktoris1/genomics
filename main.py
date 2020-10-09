from re import findall
from codes import codes, people
import json
import scipy.stats

with open('result.txt', 'r') as file:
 data = file.readlines()[1:]


data_sources = {int(pos): {int(ind): {} for ind in people} for pos in codes}
field_names = [
    'adenine',
    'cytosine',
    'guanine',
    'thymine',
    'share'
]

for i in range(len(data)//8):
    user, position = findall('[0-9]{1,50}', data[i*8])
    length = len(data_sources[int(position)][int(user)])
    for j in range(4):
       # data_sources[int(position)][int(user)].append(findall('[0-9]{1,10}', data[i*8+j+1])[0])
        data_sources[int(position)][int(user)][field_names[j]] = int(findall('[0-9]{1}.?[0-9]{0,20}', data[i*8+j+1])[0])
    data_sources[int(position)][int(user)][field_names[4]] = float(findall('[0-9]{1}.?[0-9]{0,20}', data[i * 8 + 6])[0])

    # X.append((user,position))


# with open('repeated.txt', 'w') as f:
#     for pair in set(X):
#         if X.count(pair) > 1:
#             print("number", pair[0], "in position", pair[1], "(", X.count(pair), "times )", file=f)


# with open('grouped_data.json','w') as file:
#      json.dump(data_sources,file)
number = 0
err_rate = 0.004771907590174684
rejected = 0
for position, lst in data_sources.items():
    for person, field in lst.items():
        if field.get("share", 0) == 1:
            number += 1
            GACT = [field.get('guanine'), field.get('adenine'), field.get('cytosine'), field.get('thymine')]
            n = 0
            for elem in GACT:
                n += elem
            sample_err_rate = (n - max(GACT)) / n
            diff = abs(sample_err_rate - err_rate)
            not_p_value = 0
            l = int(n * (diff - err_rate)) + 1
            r = int(n * (diff + err_rate))
            for k in range(l, r):
                not_p_value += scipy.stats.binom.pmf(k, n, err_rate)
            p_value = 1 - not_p_value
            if p_value < 0.05:
                rejected += 1

print((rejected / number ) * 100)





