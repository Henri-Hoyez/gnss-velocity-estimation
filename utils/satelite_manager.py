



def get_satelites(sats, time, reciver_informations):
    toe_limit = time + 2*(60**2)
    actual_sight = list(reciver_informations.index)
    filtered_sats = list( filter( lambda s : s.toe > time and s.toe < toe_limit and s.name[0] == 'G' and s.name in actual_sight, sats))
    
    return filtered_sats[:3]
