wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
rgi load --card_json card.json --local

