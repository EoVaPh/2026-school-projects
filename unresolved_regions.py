import glob

record_file_paths = glob.glob('extracted_chains/*')

for path in record_file_paths:
    record_file = open(path, 'r')
    records = record_file.readlines()
    record_file.close()

    unresolved_regions_count = 0

    interrupting_resnums = []

    resnum_ = None
    for record in records:
        resnum = int(record.split()[4])

        if resnum_ is not None and resnum != resnum_ + 1:
            print(path, resnum)
            interrupting_resnums.append(resnum)
            unresolved_regions_count += 1

        resnum_ = resnum

    if unresolved_regions_count > 1:
        print('ALERT', path, unresolved_regions_count, interrupting_resnums)
