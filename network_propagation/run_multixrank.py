import sys
import multixrank

multixrank_obj = multixrank.Multixrank(config=f"config.yml", wdir=f".")

ranking_df = multixrank_obj.random_walk_rank()

multixrank_obj.write_ranking(ranking_df, path=f"output")