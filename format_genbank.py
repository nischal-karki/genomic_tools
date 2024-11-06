from datetime import date

# To Do: multiple references and authors
def makegb(**kwargs) -> str:
    """
    Makes a genbank file from the keyword arguments.
    Note this should only be called by other functions or part of a python script.
    Args:
        locus: The locus of the sequence
        seq: The sequence
        seq_type: The type of sequence
        topology: The topology of the sequence
        division: The division of the sequence
        date: The date of the sequence generation
        definition: The definition of the sequence
        accession: The accession number of the sequence
        version: The version of the sequence
        keywords: The keywords of the sequence
        source: The source of the sequence
        organism: The organism of the sequence
        taxonomy: The taxonomy of the sequence
        reference: The reference of the sequence.
        authors: The authors of the sequence
        title: The title of the article
        journal: The journal of the article
        features: List of features in the format: 
            [
                {
                    "type":type, 
                    "location":location, 
                    "qualifier1":value1, 
                    "qualifier2":value2, 
                    ...
                },
            ]
    Returns:
        The genbank file as a string
    """
    kwargs = { k.lower():v for k,v in kwargs.items()}
    features = kwargs.get("features",[])
    gb = f"""LOCUS       {kwargs.get("locus","Unknown")}                {len(kwargs.get("seq",""))} {"bp" if kwargs.get("seq_type","").lower() != "protein" else "aa"} {kwargs.get("seq_type","")} {kwargs.get("topology","")} {kwargs.get("division","")} {kwargs.get("date",date.today().strftime("%d-%b-%Y"))}
DEFINITION  {kwargs.get("definition","")}
ACCESSION   {kwargs.get("accession","")}
VERSION     {kwargs.get("version","")}
KEYWORDS    {kwargs.get("keywords","")}
SOURCE      {kwargs.get("source","")}
  ORGANISM  {kwargs.get("organism","")}
            {kwargs.get("taxonomy","")}
""" + f"""REFERENCE   {kwargs.get("reference","")}
    AUTHORS   {kwargs.get("authors","")}
    TITLE     {kwargs.get("title","")}
    JOURNAL   {kwargs.get("journal","")}
""" + f"""FEATURES             LOCATION/QUALIFIERS
    source          {kwargs.get("source","")}
""" + "".join(
        [
            f'    {feature["type"]}          {feature["location"]}\n' + "\n".join(
                [f"                     /{k}={v}" for k, v in feature.items() if k != "type" and k != "location"]
            ) + "\n" for feature in features
        ]
    ) + "ORIGIN" + "\n" + "\n".join(
        [
            f"{i+1:>9}" + " " + " ".join(
                [
                    kwargs.get("seq","")[j:j+10]  
                        for j in range(i,i+60,10)
                ]
                ) for i in range(0,len(kwargs.get("seq","")),60)
        ]
    )
    return gb

def help() -> None:
    print(makegb.__doc__)

if __name__ == "__main__":
    help()