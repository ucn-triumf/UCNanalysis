# experiment-general constants used in various calculations
# Derek Fujimoto
# Aug 2024

# The raw beam on and beam off times in EPICS are in units of beam buckets, which depend on the cyclotron frequency.  We kick a fraction of the beam buckets to the UCN target; there is a blanking period between buckets when there is no beam; we use that blanking period to ramp up/down the UCN kicker magnet.  So the fundamental unit of UCN beam kicking is beam buckets.
beam_bucket_duration_s = 0.00088801 # seconds

