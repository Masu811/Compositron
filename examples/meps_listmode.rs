use compositron::{elistcol_match, eventlist::eventlist::{EventListColumnType}, importers::eventlist::bin_importer::import_bin_eventlist_row_major};

fn main() -> anyhow::Result<()> {
    let columns: Vec<(String, EventListColumnType)> = vec![
        ("e".into(), EventListColumnType::F64),
        ("t1".into(), EventListColumnType::F64),
        ("t2".into(), EventListColumnType::F64),
        ("t3".into(), EventListColumnType::F64),
    ];

    let keys = columns.iter().map(|x| x.0.clone()).collect::<Vec<_>>();

    let elist = import_bin_eventlist_row_major(
        "/mnt/YSZ-Ref-02_bin/10_YSZ-Ref-02_4.50keV_295.15K_10.0Mio_20.06.25_19.53.07-ListMode_Ch0.bin",
        columns,
        "A".into(),
        compositron::core::utils::ByteOrder::BigEndian,
        0
    )?;

    for i in 0..20 {
        for key in &keys {
            let col = elist.columns.get(key).unwrap();
            elistcol_match!(
                col,
                arr => print!("{:25} ", arr[i])
            );
        }
        println!();
    }

    Ok(())
}
